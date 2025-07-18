import sys
import os

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python get_feature_name.py <feature_file_name>")
        print("We use several external datasets in our analysis. In this script, we extract the feature "
              "represented by the file name. The script involves several hardcoded names so is not a general purpose "
              "script and is probably being used in a specific context.")
        sys.exit(1)

    # extract the filename from the path
    path = sys.argv[1]
    filename = os.path.basename(path)

    # prepare output variable
    output_feature_list = []

    if filename.startswith("GSM58566"):
        # this is probably from the aiello et al. study.
        for feature in ["Rpb1", "Rnh1", "Rnh201"]:
            if "_" + feature + "_WT_" in filename:
                output_feature_list.append(feature + "WTCRAC")

        for time_point in ["Asyn", "S30", "S45"]:
            if "_" + time_point + "_" in filename:
                output_feature_list.append(time_point)

        for replicate in ["R1", "R2"]:
            if "_" + replicate + "_" in filename:
                output_feature_list.append(replicate)

        for site in ["TES", "TSS"]:
            if site in filename:
                output_feature_list.append(site)

        if "no_zero_values" in filename:
            output_feature_list.append("no_zeroes")
        elif "zero_values" in filename:
            output_feature_list.append("zeroes")

    if len(output_feature_list) == 0:
        raise ValueError("Could not extract feature from filename: " + filename)
    else:
        print("_".join(output_feature_list))
        sys.exit(0)
