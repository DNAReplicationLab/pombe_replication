# ----------------------------------------------------------
# Copyright 2020-2021 Earlham Institute
# Written by Conrad Nieduszynski (conrad.nieduszynski@earlham.ac.uk)
#            and Sathish Thiyagarajan (sathish.thiyagarajan@earlham.ac.uk)
#            with the aid of ChatGPT and GitHub Copilot.
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
# ----------------------------------------------------------

# Usage: < path.to.detect.file python rDNA_detectSummary.py outputDirectory [n_sd_prefactor]
# example commands where the input is a mod bam file:
#   samtools view -h test_output.mod.bam | python convert_modBAM_to_detect.py | python rDNA_detectSummary.py output_dir
# example commands where the input is a detect file:
#   < test_output.detect python rDNA_detectSummary.py output_dir
#   cat test_output.detect | python rDNA_detectSummary.py output_dir

import sys
import os
import math
import re
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
from DNAscentTools.modBAM_tools import convert_detect_into_detect_stream

# check for correct number of arguments
if len(sys.argv) < 2:
    print("Usage: < path.to.detect.file python rDNA_detectSummary.py outputDirectory [n_sd_prefactor]")
    print("Example commands where the input is a mod bam file:")
    print("  samtools view -h test_output.mod.bam | python convert_modBAM_to_detect.py | python rDNA_detectSummary.py "
          "output_dir")
    print("Example commands where the input is a detect file:")
    print("  < test_output.detect python rDNA_detectSummary.py output_dir")
    print("  cat test_output.detect | python rDNA_detectSummary.py output_dir 2")
    print("Output: various files are sent to output_dir/")
    print("Optional argument: n_sd_prefactor is the factor to the s.d. by which we filter steps. "
          "Default is 3. To specify a different value, provide it as the second argument and remove "
          "the square brackets. Must be an integer (not a floating point number) because of the way we name files.")
    print("Please refer to the comments in the script for more information")
    exit(1)

# ensure that data is piped in
if sys.stdin.isatty():
    raise NotImplementedError('Please pipe in inputs')

# initialise output directory
outputDir = sys.argv[1]

# get the pre-factor for the standard deviation
if len(sys.argv) > 2:
    prefactor = int(sys.argv[2])
else:
    prefactor = 3

# create output directory if it doesn't exist
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

# set various thresholds
detectThres = 0.5  # threshold to mark T at a site as BrdU
readThres = 0.05  # threshold to mark a read as nascent
win = 290  # window size to calculate a sliding average per read
sWin = win * 3  # window size used while calculating BrdU steps between two consecutive windows
lWin = sWin * 5  # two detected peaks in step size vs coordinate should be at least this distance apart.
# set to 1 to disable as that is the lowest possible distance between two peaks anyway.
# theoretically, based on our filter structure, this should be set to some factor times sWin,
# depending on our tolerance.

step_size_thres_3sigma_th = prefactor * math.sqrt(2) * 0.5 / math.sqrt(sWin)
# theoretical threshold for step size, read below
# we can ask: let's say we have pure, uncorrelated, random data i.e. equal probability of zero and one at
# each thymidine with no memory from one thymidine to the next. Then, what do we expect the step size to be?
# Using stats, we can calculate the standard deviation in step size is sqrt(2) * 0.5 / sqrt(window size).
# We ask that any step we detect is at least N standard deviations significant w.r.t. a random process.

# initialise various buffers
stepBuff = []
stepCountBuff = []
stepBEDbuff = []
nascBEDbuff = []
parBEDbuff = []
posBuff = []
leftPosBuff = []
rightPosBuff = []
leftDensityBuff = []
rightDensityBuff = []
detectIndexBuff = []

# initialise various counters
readCount = 0
nascentCount = 0

# initialise header line
txtHeader = '\t'.join(['chr',
                       'start',
                       'end',
                       'read',
                       '[BrdU]',
                       'strand',
                       'mapped length',
                       'highest BBT',
                       'largest step',
                       'start-end diff'])


def processRead(detectThreshold: float, readThreshold: float, Window: int, sWindow: int, lWindow: int,
                chromosome: str, readID: str, strand: str, mappingStart: int,
                mappingEnd: int, coordBuff: list[int], detectBuff: list[float],
                stepSizeThres: float) \
        -> tuple[str, list, list, tuple, list, list, list, list, list, str]:
    r""" Process a read and return step information and whether it is a nascent read

    Args:
        detectThreshold: get hard BrdU calls (0 or 1) from soft calls (probability b/w 0 and 1) using this threshold
        readThreshold: read marked as nascent if the average BrdU content is above this threshold
        Window: window size to calculate a sliding average per read
        sWindow: window size to calculate BrdU steps between two consecutive windows
        lWindow: peaks in step size vs coordinate should be at least this distance apart
        chromosome: chromosome/contig of the read
        readID: read ID of the read
        strand: "+" or "-" depending on the direction of the alignment of the read to the reference
        mappingStart: start position of the alignment on the reference
        mappingEnd: end position of the alignment on the reference
        coordBuff: list of coordinates of the BrdU soft calls
        detectBuff: list of BrdU soft calls
        stepSizeThres: peaks in step size below this threshold are rejected.

    Returns:
        A tuple containing the following:
            - 'NA'/'Nascent'/'Parent' depending on whether the read is not long enough, is a nascent or a parental read
            - []/[]/list of all step sizes in the read corresponding to the three cases above
            - []/[]/list of local peaks in step size in 6+1 bed format corresponding to the three cases above
            - ()/tuple of 10 pieces of data in 6+4 bed format/tuple of 12 pieces of data in 6+4+2 bed format
              of read level data corresponding to the three cases above
            - []/[]/list of coordinates corresponding to the step sizes in item 2 above
            - []/[]/list of left-side values corresponding to the step sizes in item 2 above
            - []/[]/list of right-side values corresponding to the step sizes in item 2 above
            - []/[]/list of coordinates corresponding to the left most point of step size calculation in item 2 above
            - []/[]/list of coordinates corresponding to the right most point of step size calculation in item 2 above
            - detectIndex: a unique read identifier in the format readID_chromosome_mappingStart_mappingEnd_strand

    Note:
        BED format generally means six pieces of data which will be separated by tabs when they are output in a file.
        The six pieces of data are: chromosome, start, end, name, score, strand.
        Score is a number between 0 and 1000, but we use other numbers here.
        The 6+N bed format is the same as the 6 bed format, but with N additional columns of data.
        The 6+4+2 bed format is something we made up and is not a standard format.

    """
    # construct the detectIndex
    detect_index = f"{readID}_{chromosome}_{mappingStart}_{mappingEnd}_" + ("fwd" if strand == "+" else "rev")

    # Create a numpy array from the list `detectBuff` and store it in `BrdU_array`
    BrdU_array = np.array(detectBuff)

    # if the data is too short, skip it
    if BrdU_array.size <= (2 * sWindow):
        return 'NA', [], [], (), [], [], [], [], [], detect_index

    # Create a numpy array from the list `coordBuff` and store it in `coordArray`
    coordArray = np.array(coordBuff)

    # get hard BrdU calls by thresholding
    hardB = np.array(BrdU_array >= detectThreshold, dtype=float)

    # Calculate the difference between two uniform filtered arrays of `BrdU_array` with sizes `sWindow`
    left_average = uniform_filter1d(hardB, size=sWindow, mode='wrap', origin=(sWindow // 2) - 1)
    right_average = uniform_filter1d(hardB, size=sWindow, mode='wrap', origin=-((sWindow // 2) - 1))
    sWindowDiff = left_average - right_average

    # Get the absolute value of the `sWindowDiff`
    absWindowDiff = abs(sWindowDiff)

    # calculate mean brdu for the entire read
    readBBT = np.sum(hardB) / hardB.size

    # Filter the `hardB` array using a uniform filter of size `Window`
    BBT = uniform_filter1d(hardB, size=Window, origin=0)

    # create a temporary output object
    BEDtemp = (chromosome,
               mappingStart,
               mappingEnd,
               readID,
               readBBT,
               strand,
               (mappingEnd - mappingStart),
               BBT[Window:-Window].max(),
               absWindowDiff[sWindow:-sWindow].max(),
               absWindowDiff[len(absWindowDiff) - 1])

    # Check if the read is a nascent read
    if readBBT >= readThreshold:
        # Find the indices of peaks in `absWindowDiff[sWindow:-sWindow]` with minimum distance and minimum height
        peakInStepIndices, _ = find_peaks(absWindowDiff[sWindow:-sWindow], distance=lWindow, height=stepSizeThres)

        # Get the corresponding coordinate values from `coordArray`
        peakInStepCoords = coordArray[peakInStepIndices + sWindow]

        # Get the step size values from `absWindowDiff`
        peakInStepSizes = absWindowDiff[peakInStepIndices + sWindow]

        stepBEDList = [(chromosome, peakInStepCoords[i], peakInStepCoords[i],
                        readID, float(peakInStepSizes[i]), strand, (mappingEnd - mappingStart))
                       for i in range(len(peakInStepSizes))]

        return 'Nascent', sWindowDiff[sWindow:-sWindow].tolist(), stepBEDList, \
            (*BEDtemp, peakInStepCoords, peakInStepSizes), coordArray[sWindow:-sWindow].tolist(), \
            left_average[sWindow:-sWindow].tolist(), right_average[sWindow:-sWindow].tolist(), \
            coordArray[:-2 * sWindow].tolist(), coordArray[2 * sWindow:].tolist(), detect_index

    else:
        return 'Parent', [], [], BEDtemp, [], [], [], [], [], detect_index


def convert_lists_to_strings(np_nd_array: np.ndarray) -> np.ndarray:
    """ Convert all 1D arrays in a 2D numpy array into strings.

    Args:
        np_nd_array: A 2D numpy array.

    Returns:
        A 2D numpy array with all 1D arrays converted into strings.
    """

    # Check if the input is a 2D array
    if np_nd_array.ndim != 2:
        raise ValueError("The input should be a 2D array.")

    # Loop through each element of the np_nd_array
    for i in range(np_nd_array.shape[0]):
        for j in range(np_nd_array.shape[1]):

            # Check if the element is of type ndarray and dimension 1
            if isinstance(np_nd_array[i, j], np.ndarray) and np_nd_array[i, j].ndim == 1:
                # Convert the 1D array into a string as specified
                str_representation = "[" + " ".join(str(elem) for elem in np_nd_array[i, j]) + "]"
                np_nd_array[i, j] = re.sub(r'(-?\d+\.\d+)', lambda x: '{0:.{1}f}'.format(float(x.group(1)), 4),
                                           str_representation)

    return np_nd_array


# process detect data
for k in convert_detect_into_detect_stream(sys.stdin):

    if "comments" in k:
        # skip the header of the detect file which has comments
        continue
    else:

        # update counters
        readCount += 1

        # process the read
        read_type, all_step_sizes, \
            critical_steps_bed_data, \
            read_level_bed_data, \
            pos_data, \
            left_density_data, \
            right_density_data, \
            left_position_data, \
            right_position_data, \
            index \
            = processRead(detectThres, readThres, win, sWin, lWin,
                          k["refContig"], k["readID"], "+" if k["strand"] == "fwd" else "-",
                          k["refStart"], k["refEnd"], k["posOnRef"], k["probBrdU"], step_size_thres_3sigma_th)

        # store the read in the appropriate buffer depending on its type
        if read_type == 'Nascent':
            stepBuff.extend(all_step_sizes)
            posBuff.extend(pos_data)
            stepCountBuff.append(len(all_step_sizes))
            stepBEDbuff.extend(critical_steps_bed_data)
            nascBEDbuff.append(read_level_bed_data)
            leftDensityBuff.extend(left_density_data)
            rightDensityBuff.extend(right_density_data)
            leftPosBuff.extend(left_position_data)
            rightPosBuff.extend(right_position_data)
            detectIndexBuff.append(index)
        elif read_type == 'Parent':
            parBEDbuff.append(read_level_bed_data)

stepArray = np.array(stepBuff)
stepSD = stepArray.std()
stepMean = stepArray.mean()

# If there are any steps for a nascent read, put the maximum step size in `maxStepArray`
maxStepArray = np.array([max(k[11], default=0) for k in nascBEDbuff])

# Put the absolute value of the start end analogue density difference in `seStepArray`
seStepArray = np.array([k[9] for k in nascBEDbuff])

# create various arrays and normalize step sizes by the step SD
stepBEDarray = np.array(stepBEDbuff, dtype=object)
stepBEDarray[:, 4] = stepBEDarray[:, 4] / stepSD

nascBEDarray = np.array(nascBEDbuff, dtype=object)
nascBEDarray[:, 11] = nascBEDarray[:, 11] / stepSD

# collect all step sizes of nascent reads
nascStepArray = np.concatenate(nascBEDarray[:, 11], axis=None)

# get an array of the data from the reads identified to be parental
parBEDarray = np.array(parBEDbuff, dtype=object)

print(nascStepArray[:10])

# print some statistics
print('## Total reads:', readCount)
print('## Nascent reads:', len(nascBEDbuff), ' ( Total BrdU fraction in read >', readThres, ')')
print('## Total length of nascent reads:', nascBEDarray[:, 6].sum(axis=0))
print('##')
print('## SD in steps:', str("%.3f" % stepSD), ' (', sWin, 'T positions )')
print('## Mean in steps:', str("%.3f" % stepMean), ' (', sWin, 'T positions )')
print('## (abs) Max step size:', str("%.3f" % (abs(stepArray).max())), ' (', sWin, 'T positions )')
print('## (abs) Min step size:', str("%.3f" % (abs(stepArray).min())), ' (', sWin, 'T positions )')
print('##')
# These calculations are relative to zero, but they should be relative to the mean
# The mean is very close to zero, so probably won't make much difference
# To correct this I'll need to use the step values before abs
print('## Reads with step > 1 sd:', np.sum(maxStepArray >= stepSD))
print('## Reads with step > 2 sd:', np.sum(maxStepArray >= 2 * stepSD))
print('## Reads with step > 3 sd:', np.sum(maxStepArray >= 3 * stepSD))
print('## Reads with step > 4 sd:', np.sum(maxStepArray >= 4 * stepSD))
print('##')
print('## Number of steps > 2 sd:', np.sum(nascStepArray >= 2))
print('## Number of steps > 3 sd:', np.sum(nascStepArray >= 3))
print('## Number of steps > 4 sd:', np.sum(nascStepArray >= 4))
print('##')
print('## Reads with start-end diff > 1 sd:', np.sum(seStepArray >= stepSD))
print('## Reads with start-end diff > 2 sd:', np.sum(seStepArray >= 2 * stepSD))
print('## Reads with start-end diff > 3 sd:', np.sum(seStepArray >= 3 * stepSD))
print('## Reads with start-end diff > 4 sd:', np.sum(seStepArray >= 4 * stepSD))

# identify outliers in step size
boolSEstepOutlier = seStepArray > prefactor * stepSD
boolMaxStepOutlier = maxStepArray > prefactor * stepSD

# set various formats for saving data
bedFmtString = '%s\t%d\t%d\t%s\t%.3f\t%s'
bed6Plus1FmtString = bedFmtString + '\t%d'
bed6Plus4FmtString = bed6Plus1FmtString + '\t%.3f\t%.3f\t%.3f'
bed6Plus4PlusListsFmtString = bed6Plus4FmtString + '\t%s\t%s'
bed6Plus5FmtString = bedFmtString + '\t%.3f\t%.3f\t%d\t%d\t%s'

# save text and bed format of nascent data
if nascBEDarray.tolist():
    np.savetxt(f'{outputDir}/nascent.txt', convert_lists_to_strings(nascBEDarray), fmt=bed6Plus4PlusListsFmtString,
               header=txtHeader)
    np.savetxt(f'{outputDir}/nascent.bed', nascBEDarray[:, :6], fmt=bedFmtString)
else:
    raise ImplementationError("No nascent reads detected; don't know what to do. Exiting.")

# save text and bed format of parental data
if parBEDarray.tolist():
    np.savetxt(f'{outputDir}/parental.txt', parBEDarray, fmt=bed6Plus4FmtString, header=txtHeader)
    np.savetxt(f'{outputDir}/parental.bed', parBEDarray[:, :6], fmt=bedFmtString)
else:
    print("No parental reads detected.")

# save text and bed format of start end difference data that lie outside the N SD range
np.savetxt(f'{outputDir}/start_end_diff.txt', convert_lists_to_strings(nascBEDarray[boolSEstepOutlier]),
           fmt=bed6Plus4PlusListsFmtString, header=txtHeader)
np.savetxt(f'{outputDir}/start_end_diff.bed', nascBEDarray[boolSEstepOutlier, :6], fmt=bedFmtString)

# save text and bed format of step data that lie outside the N SD range
np.savetxt(f'{outputDir}/{prefactor}SD_step.txt', convert_lists_to_strings(nascBEDarray[boolMaxStepOutlier]),
           fmt=bed6Plus4PlusListsFmtString, header=txtHeader)
np.savetxt(f'{outputDir}/{prefactor}SD_step.bed', nascBEDarray[boolMaxStepOutlier, :6], fmt=bedFmtString)

# save list of all peaks in step size per nascent read in bed format
np.savetxt(f'{outputDir}/nascentStepList.txt', stepBEDarray, fmt=bed6Plus1FmtString)

# save list of all step sizes per nascent read in bed format.
stepArrayIter = iter(stepArray)
posBuffIter = iter(posBuff)
leftPosBuffIter = iter(leftPosBuff)
rightPosBuffIter = iter(rightPosBuff)
leftDensityBuffIter = iter(leftDensityBuff)
rightDensityBuffIter = iter(rightDensityBuff)
with open(f'{outputDir}/nascentAllStepListRaw.bed', 'w') as f:
    for stepCount, nascBED, currentDetectIndex in zip(stepCountBuff, nascBEDbuff, detectIndexBuff):
        for _, currentStep, currentStart, \
                currentLeftDensity, currentRightDensity, \
                currentLeftBoundary, currentRightBoundary in zip(range(stepCount), stepArrayIter,
                                                                 posBuffIter, leftDensityBuffIter,
                                                                 rightDensityBuffIter, leftPosBuffIter,
                                                                 rightPosBuffIter):
            f.write(bed6Plus5FmtString % (nascBED[0], currentStart, currentStart + 1, nascBED[3],
                                          currentStep, nascBED[5], currentLeftDensity,
                                          currentRightDensity, currentLeftBoundary,
                                          currentRightBoundary, currentDetectIndex) + '\n')
