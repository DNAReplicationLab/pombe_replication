#!/usr/bin/env bash
#SBATCH --mem=200G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J dnascentQC
#SBATCH --mail-type=END,FAIL
#SBATCH --constraint=""
#SBATCH --time=23:59:59

# goal of program
# ================
# Perform analysis using read lengths, fork lengths, analogue density etc.
# from sequencing summary files, mod bam files, forksense files etc. and produce output in a pdf file
# (some of these inputs may be optional, read the usage section below).
# The goal is to quickly compute statistics and see if the analogue-incorporation experiment is satisfactory.
# Any sophisticated or detailed analysis will not be performed here.

# program usage
# =============
# all of the following usage examples are valid
# sbatch dnascent_qc.sh sequencing_summary mod_bam forksense_dir op_dir tag
# sbatch dnascent_qc.sh sequencing_summary mod_bam forksense_dir op_dir
# sbatch dnascent_qc.sh sequencing_summary mod_bam "" op_dir
# sbatch dnascent_qc.sh "" mod_bam "" op_dir
# sbatch dnascent_qc.sh "" mod_bam forksense_dir op_dir
# can use bash in place of sbatch
# sequencing_summary (optional): the sequencing summary file produced by the basecaller.
#                                if you don't have this file, use an empty string or pass inputs like /dev/null
# mod_bam: bam file with modification probabilities. The pipeline produces this.
# forksense_dir (optional): directory that contains files with left forks, right forks, origins, terminations
#                           information. These files are produced by forkSense in the pipeline.
#                           If you do not have this information or want to avoid forksense-related statistics,
#                           pass an empty string or some other invalid input like /dev/null.
# op_dir: output directory where the qc pdf file goes. A few other files are produced as well.
# tag (optional): thymidine modification tag. if not given or "", defaults to 'T'. do not pass invalid inputs
#                 like /dev/null here.

# load variables from the command line
sequencing_summary=
if [[ -f "$1" ]]; then
  sequencing_summary="$1";
fi
mod_bam=$2;
forksense_dir=
if [[ -d "$3" ]]; then
  forksense_dir=$(cd "$3" || exit; pwd);
fi
mkdir -p "$4";
op_dir=$(cd "$4" || exit; pwd) # set output directory
tag=${5:-T}

# set output tex file
# the output pdf is the output tex file, but with '.tex' replaced with '.pdf'
op_tex="$op_dir"/dnascent_qc.tex

# make directory where output of each step goes
analysis_op_dir="$op_dir"/dnascent_qc_data
mkdir -p "$analysis_op_dir"

# load packages
source load_package.sh -miller -latex -python -samtools -R -bedtools -modkit;

# load configuration variables
source config.sh

# create a word count function
wc_file()
{
  cat | wc -l
}

# create functions for commonly used latex tasks
page_break()
{
    echo "\pagebreak"
}

section()
{
  echo "\section{$1}"
}

sub_section()
{
  echo "\subsection{$1}"
}

begin_verbatim()
{
  echo "\begin{verbatim}"
}

end_verbatim()
{
  echo "\end{verbatim}"
}

figure()
{
  echo "\begin{figure}[h!]"
  echo "\centering"
  echo "\includegraphics[scale=0.5]{$1}"
  echo "\end{figure}"
}

# write the latex output
{

  # make initial few pages like title, table of contents etc.
  echo "\documentclass{article}"
  echo "\usepackage{graphicx}"
  echo "\usepackage{hyperref}"
  echo "\usepackage[a4paper, margin=1in]{geometry}"
  echo "\title{QC on DNAscent and other analyses from nanopore data}"
  echo "\\author{${config[name]} \\\\ ${config[email]}}"
  echo "\date{\today}"
  echo "\begin{document}"
  echo "\maketitle"
  echo "\tableofcontents"
  page_break

  # list source files used for analysis
  section "Source files/directories"
  begin_verbatim
  {
    echo "Sequencing summary: $sequencing_summary";
    echo "Mod bam file: $mod_bam";
    echo "Forksense directory: $forksense_dir";
  } | fold -w 80
  end_verbatim
  page_break

  if [ ! "$sequencing_summary" == "" ]; then
    section "Read lengths"
    sub_section "Data for histogram of read lengths from sequencing summary file"
    begin_verbatim

    # raw data for histogram of read lengths
    # shellcheck disable=SC1010
    mlr --itsv --opprint cut -f sequence_length_template then histogram -f sequence_length_template\
      --lo 0 --hi 200000 --nbins 40 "$sequencing_summary" | tee "$analysis_op_dir"/histogram_read_lengths

    end_verbatim
    page_break

    label=histogram_read_lengths
    < "$analysis_op_dir"/"$label" \
        Rscript plotting_and_short_analyses/plot_histogram.R sequence_length_template_count \
           "$analysis_op_dir"/"$label"_plot.png "Read length (b)"

    figure "$analysis_op_dir"/"$label"_plot.png
    page_break

    sub_section "Data for yield in bases vs binned read length from sequencing summary file"
    begin_verbatim

    # there are many occurrences of a number like 5000 in the command below
    # that's the bin size. if you want some other bin size, change all instances of this number
    # in the command below.
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    mlr --itsv --opprint rename sequence_length_template,l \
      then cut -f l then put '$p=round($l/5000)' \
      then stats1 -a sum -f l -g p then put '$bin_lo = $p * 5000; $bin_hi = ($p + 1) * 5000' \
      then sort -n p then cut -o -f bin_lo,bin_hi,l_sum then rename l_sum,num_bases "$sequencing_summary" |\
      tee "$analysis_op_dir"/histogram_yield_data

    end_verbatim
    page_break

    label=histogram_yield_data
    < "$analysis_op_dir"/"$label" \
        Rscript plotting_and_short_analyses/plot_histogram.R num_bases \
           "$analysis_op_dir"/"$label"_plot.png "Read length (b)" "Yield (b)"

    figure "$analysis_op_dir"/"$label"_plot.png
    page_break

    sub_section "Statistics of read lengths from sequencing summary file"
    begin_verbatim

    # read length statistics
    mlr --itsv --oxtab stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f sequence_length_template\
      "$sequencing_summary"
    echo -e "\n N50 (sampling 100,000 reads randomly) \n"
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    mlr --itsv --oxtab shuffle then head -n 100000 then cut -f sequence_length_template then \
        sort -n sequence_length_template \
        then fraction -f sequence_length_template -c then filter '$sequence_length_template_cumulative_fraction > 0.5' \
        then head -n 1 "$sequencing_summary"

    end_verbatim
    page_break

    sub_section "Read length from sequencing summary file vs alignment length from mod bam file"
    sbatch --wait join_alignLen_basecalledLen_per_readid.sh "$sequencing_summary" "$mod_bam" \
      "$analysis_op_dir"/table_basecalled_align_len_read_id > /dev/null
    # we wait here till the job finishes.
    # shellcheck disable=SC2016
    < "$analysis_op_dir"/table_basecalled_align_len_read_id \
        mlr --itsv --otsv --skip-comments put '$align_length = $align_length/1000;
            $sequence_length_template = $sequence_length_template/1000' |\
        Rscript plotting_and_short_analyses/plot_scatter.R align_length sequence_length_template \
           "$analysis_op_dir"/plot_table_basecalled_align_len_read_id.png "Alignment length (kb)" \
            "Basecalled length (kb)" 0,100 0,100 1 0
    figure "$analysis_op_dir"/plot_table_basecalled_align_len_read_id.png
    begin_verbatim
    echo "NOTE: The plot runs from x = 0 to x = 100 and same for y. Data outside this range are not shown."
    end_verbatim
    page_break
  fi

  if [ ! "$forksense_dir" == "" ]; then
    section "Fork and origin statistics"
    sub_section "Numbers of different features"
    begin_verbatim

    echo "Number of left forks"
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed wc_file
    echo "Number of right forks"
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed wc_file
    echo "Number of origins"
    < "$forksense_dir"/origins_DNAscent_forkSense.bed wc_file
    echo "Number of terminations"
    < "$forksense_dir"/terminations_DNAscent_forkSense.bed wc_file
    echo "Number of molecules with left forks"
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed awk '{print $4}' | sort | uniq | wc -l
    echo "Number of molecules with right forks"
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed awk '{print $4}' | sort | uniq | wc -l

    echo ""
    echo "Considering only fwd reads"
    echo "Number of left forks"
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed grep fwd | wc_file
    echo "Number of right forks"
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed grep fwd | wc_file
    echo "Number of origins"
    < "$forksense_dir"/origins_DNAscent_forkSense.bed grep fwd | wc_file
    echo "Number of terminations"
    < "$forksense_dir"/terminations_DNAscent_forkSense.bed grep fwd | wc_file
    echo "Number of molecules with left forks"
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed grep fwd | awk '{print $4}' | sort | uniq | wc -l
    echo "Number of molecules with right forks"
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed grep fwd | awk '{print $4}' | sort | uniq | wc -l

    echo ""
    echo "Considering only rev reads"
    echo "Number of left forks"
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed grep rev | wc_file
    echo "Number of right forks"
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed grep rev | wc_file
    echo "Number of origins"
    < "$forksense_dir"/origins_DNAscent_forkSense.bed grep rev | wc_file
    echo "Number of terminations"
    < "$forksense_dir"/terminations_DNAscent_forkSense.bed grep rev | wc_file
    echo "Number of molecules with left forks"
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed grep rev | awk '{print $4}' | sort | uniq | wc -l
    echo "Number of molecules with right forks"
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed grep rev | awk '{print $4}' | sort | uniq | wc -l

    end_verbatim
    page_break

    sub_section "Raw data for histogram of fork lengths"
    begin_verbatim

    # raw data for histogram of fork lengths
    echo -e "\n All forks\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    cat "$forksense_dir"/leftForks_DNAscent_forkSense.bed "$forksense_dir"/rightForks_DNAscent_forkSense.bed |\
    mlr --ipprint --opprint --implicit-csv-header put '$fork_length = $3-$2' then histogram -f fork_length \
      --lo 0 --hi 200000 --nbins 40 | tee "$analysis_op_dir"/histogram_all_fork_lengths

    echo -e "\n Left forks\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed \
    mlr --ipprint --opprint --implicit-csv-header put '$fork_length = $3-$2' then histogram -f fork_length \
      --lo 0 --hi 200000 --nbins 40 | tee "$analysis_op_dir"/histogram_left_fork_lengths

    echo -e "\n Right forks\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed \
    mlr --ipprint --opprint --implicit-csv-header put '$fork_length = $3-$2' then histogram -f fork_length \
      --lo 0 --hi 200000 --nbins 40 | tee "$analysis_op_dir"/histogram_right_fork_lengths

    end_verbatim
    page_break

    label=histogram_all_fork_lengths
    < "$analysis_op_dir"/"$label" \
        Rscript plotting_and_short_analyses/plot_histogram.R fork_length_count \
           "$analysis_op_dir"/"$label"_plot.png "All fork length (b)"

    figure "$analysis_op_dir"/"$label"_plot.png
    page_break

    label=histogram_left_fork_lengths
    < "$analysis_op_dir"/"$label" \
        Rscript plotting_and_short_analyses/plot_histogram.R fork_length_count \
           "$analysis_op_dir"/"$label"_plot.png "Left fork length (b)"

    figure "$analysis_op_dir"/"$label"_plot.png
    page_break

    label=histogram_right_fork_lengths
    < "$analysis_op_dir"/"$label" \
        Rscript plotting_and_short_analyses/plot_histogram.R fork_length_count \
           "$analysis_op_dir"/"$label"_plot.png "Right fork length (b)"

    figure "$analysis_op_dir"/"$label"_plot.png
    page_break

    sub_section "Statistics of fork lengths"
    begin_verbatim

    # fork length statistics
    echo -e "\n All forks\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    cat "$forksense_dir"/leftForks_DNAscent_forkSense.bed "$forksense_dir"/rightForks_DNAscent_forkSense.bed |\
    mlr --ipprint --oxtab --implicit-csv-header put '$fork_length = $3-$2' then \
      stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f fork_length

    echo -e "\n N50 \n"
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    cat "$forksense_dir"/leftForks_DNAscent_forkSense.bed "$forksense_dir"/rightForks_DNAscent_forkSense.bed |\
    mlr --ipprint --oxtab --implicit-csv-header put '$fork_length = $3-$2' then \
        cut -f fork_length then sort -n fork_length \
        then fraction -f fork_length -c then filter '$fork_length_cumulative_fraction > 0.5' \
        then head -n 1

    echo -e "\n Left forks\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed \
    mlr --ipprint --oxtab --implicit-csv-header put '$fork_length = $3-$2' then \
      stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f fork_length

    echo -e "\n N50 \n"
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    < "$forksense_dir"/leftForks_DNAscent_forkSense.bed \
    mlr --ipprint --oxtab --implicit-csv-header put '$fork_length = $3-$2' then \
        cut -f fork_length then sort -n fork_length \
        then fraction -f fork_length -c then filter '$fork_length_cumulative_fraction > 0.5' \
        then head -n 1

    echo -e "\n Right forks\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed \
    mlr --ipprint --oxtab --implicit-csv-header put '$fork_length = $3-$2' then \
      stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f fork_length

    echo -e "\n N50 \n"
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    < "$forksense_dir"/rightForks_DNAscent_forkSense.bed \
    mlr --ipprint --oxtab --implicit-csv-header put '$fork_length = $3-$2' then \
        cut -f fork_length then sort -n fork_length \
        then fraction -f fork_length -c then filter '$fork_length_cumulative_fraction > 0.5' \
        then head -n 1

    end_verbatim
    page_break

    sub_section "Statistics of origin, termination uncertainty intervals"
    echo "DNAscent associates one genomic window per called origin/termination."
    echo "Here are the statistics of those windows."
    begin_verbatim

    echo -e "\n All origins\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    < "$forksense_dir"/origins_DNAscent_forkSense.bed \
    mlr --ipprint --oxtab --implicit-csv-header put '$origin_length = $3-$2' then \
      stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f origin_length

    echo -e "\n All terminations\n"
    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    < "$forksense_dir"/terminations_DNAscent_forkSense.bed \
    mlr --ipprint --oxtab --implicit-csv-header put '$termination_length = $3-$2' then \
      stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f termination_length

    end_verbatim
    page_break
  fi

  section "Whole read analogue density statistics"
  sub_section "Raw data for histogram of densities"
  begin_verbatim

  # calculate mean analogue density over whole reads
  # we are appropriating a nascent read identifying script for this purpose
  samtools view "$mod_bam" |\
    python get_modBAM_nascent_reads.py --window 0 --nascentThres 0 --showVal --tag "$tag" > "$mod_bam"_statistics

  echo -e "\n All molecules\n"
  # shellcheck disable=SC1010
  # shellcheck disable=SC2016
  < "$mod_bam"_statistics \
  mlr --itsv --opprint --implicit-csv-header label read_id,mean_brdU \
      then histogram -f mean_brdU --lo 0 --hi 1 --nbins 50 |\
      tee "$analysis_op_dir"/histogram_analogue_densities

  end_verbatim
  page_break

  label=histogram_analogue_densities
  < "$analysis_op_dir"/"$label" \
      Rscript plotting_and_short_analyses/plot_histogram.R mean_brdU_count \
         "$analysis_op_dir"/"$label"_plot.png "Analogue density"

  figure "$analysis_op_dir"/"$label"_plot.png
  page_break

  sub_section "Analogue density statistics"
  begin_verbatim

  echo -e "\n All molecules\n"
  # shellcheck disable=SC1010
  # shellcheck disable=SC2016
  < "$mod_bam"_statistics \
  mlr --itsv --oxtab --implicit-csv-header label read_id,mean_brdU \
      then stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f mean_brdU

  end_verbatim
  page_break

  sub_section "Statistics of reads and read lengths in modbam file"
  begin_verbatim

  # make a temporary file
  temp_file_mod_bam_alignment_info=$(mktemp)

  # read length statistics
  # shellcheck disable=SC1010
  bedtools bamtobed -i "$mod_bam" > "$temp_file_mod_bam_alignment_info"

  # shellcheck disable=SC1010
  < "$temp_file_mod_bam_alignment_info" \
    awk '{print $3-$2}' |\
    mlr --itsv --oxtab --implicit-csv-header label l then stats1 -a count,sum,min,p10,p50,mean,p90,max,stddev -f l

  samtools stats "$mod_bam" > "$mod_bam"_statistics_2

  echo -e "\n N50 of reads within modbam file \n"
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$mod_bam"_statistics_2 grep ^RL | cut -f 2- |\
    mlr --itsv --opprint --implicit-csv-header label l,count \
        then put '$b = $l * $count' then fraction -f b -c \
        then filter '$b_cumulative_fraction > 0.5' \
        then cut -f l,b_cumulative_fraction \
        then rename b_cumulative_fraction,l_cumulative_fraction \
        then head -n 1

  echo -e "\n Number of forward and reverse reads within modbam file \n"
  # count number of forward and reverse reads
  # shellcheck disable=SC1010
  < "$temp_file_mod_bam_alignment_info" \
    mlr --tsv --implicit-csv-header cut -f 6 then rename 6,orientation then count -g orientation

  # remove temporary file
  rm "$temp_file_mod_bam_alignment_info"

  end_verbatim
  page_break

  sub_section "Raw data for yield in bases vs binned read length from mod bam file"
  begin_verbatim

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$mod_bam"_statistics_2 grep ^RL | cut -f 2- |\
    mlr --itsv --opprint --implicit-csv-header label l,count \
        then put '$b = $l * $count' \
        then put '$p=round($l/5000)' \
        then cut -f b,p \
        then stats1 -a sum -f b -g p \
        then put '$bin_lo = $p * 5000; $bin_hi = ($p + 1) * 5000' \
        then sort -n p \
        then cut -o -f bin_lo,bin_hi,b_sum \
        then rename b_sum,n_bases |\
        tee "$analysis_op_dir"/histogram_yield_data_from_mod_bam

  end_verbatim
  page_break

  label=histogram_yield_data_from_mod_bam
  < "$analysis_op_dir"/"$label" \
      Rscript plotting_and_short_analyses/plot_histogram.R n_bases \
         "$analysis_op_dir"/"$label"_plot.png "Read length (b)" "Yield (b)"

  figure "$analysis_op_dir"/"$label"_plot.png
  page_break

  section "Windowed analogue density statistics"
  sub_section "Raw data for histogram of densities"
  begin_verbatim

  echo "Using a window size of 300 thymidines."
  echo "NOTE: As this calculation is compute-intensive, we choose 5% of the reads at random and calculate the windowed analogue density along them." |\
    fold -w 80

  # calculate mean analogue density over windows in reads
  {
    echo "# windowed analogue densities of a random subset of reads in 300T windows";

    samtools view -s 0.05 "$mod_bam" |\
      python get_modBAM_windowed_data.py --window 300 --tag "$tag";

  }  > "$mod_bam"_statistics_3

  echo -e "\n Subset of molecules\n"
  # shellcheck disable=SC1010
  # shellcheck disable=SC2016
  < "$mod_bam"_statistics_3 \
  mlr --itsv --opprint --skip-comments --implicit-csv-header label read_id,mean_brdU \
      then histogram -f mean_brdU --lo 0 --hi 1 --nbins 50 |\
      tee "$analysis_op_dir"/histogram_windowed_analogue_densities

  end_verbatim
  page_break

  label=histogram_windowed_analogue_densities
  < "$analysis_op_dir"/"$label" \
      Rscript plotting_and_short_analyses/plot_histogram.R mean_brdU_count \
         "$analysis_op_dir"/"$label"_plot.png "Analogue density"

  figure "$analysis_op_dir"/"$label"_plot.png
  page_break

  section "Raw analogue probability statistics"
  sub_section "Raw data for histogram of probabilities"
  begin_verbatim

  echo "NOTE: As this calculation is compute-intensive, we choose 10000 reads at random" | fold -w 80
  echo "NOTE: This calculation may count bases whose modification status is unknown as unmodified" | fold -w 80

  echo -e "\n Subset of molecules\n"

  input_mod_bam_adjust_tag=${mod_bam}.tag_adjusted.bam
  histogram_folder=$(dirname "$(realpath "$mod_bam")")/modkit_histogram

  # convert mod BAM from using T+T to T+B
  modkit adjust-mods --convert T B "$mod_bam" "$input_mod_bam_adjust_tag"
  samtools index "$input_mod_bam_adjust_tag"

  # sample probabilities and make a histogram
  modkit sample-probs -p 0.1,0.2,0.3,0.4,0.5 \
    "$input_mod_bam_adjust_tag" -o "$histogram_folder" --hist -n 10000 --force

  rm "$input_mod_bam_adjust_tag"
  rm "$input_mod_bam_adjust_tag".bai

  probabilities_tsv="$histogram_folder"/probabilities.tsv
  temp_file_probabilities="$histogram_folder"/temp_file_probabilities
  # shellcheck disable=SC1010
  # shellcheck disable=SC2016
  mlr --tsv put \
    'if($code == "-" && $primary_base == "T"){
      $bin_lo = 1 - $range_end;
      $bin_hi = 1 - $range_start;
    }elif($code == "B" && $primary_base == "T"){
      $bin_lo = $range_start;
      $bin_hi = $range_end;
    }else {
      # this is an error state.
      # we are not supposed to enter here.
    }'\
   then cut -f bin_lo,bin_hi,count then sort -n bin_lo "$probabilities_tsv" > "$temp_file_probabilities"

  cat "$temp_file_probabilities"

  end_verbatim
  page_break

  label=histogram_analogue_probabilities
  < "$temp_file_probabilities" \
      Rscript plotting_and_short_analyses/plot_histogram.R count \
         "$analysis_op_dir"/"$label"_plot.png "Prob modification per base" "Count" 0,1 0,auto

  figure "$analysis_op_dir"/"$label"_plot.png
  page_break

  echo -e "\n Zoom in to probabilities greater than 0.04\n"

  # shellcheck disable=SC2016
  mlr --tsv filter '$bin_lo > 0.04' "$temp_file_probabilities" |\
    Rscript plotting_and_short_analyses/plot_histogram.R count \
         "$analysis_op_dir"/"$label"_plot_zoom.png "Prob modification per base" "Count" 0,1 0,auto

  figure "$analysis_op_dir"/"$label"_plot_zoom.png
  page_break

  rm "$temp_file_probabilities"

  echo "\end{document}"

} > "$op_tex"

# convert to pdf
# pdflatex sometimes needs to be run twice, so we are going to run it twice every time. This is not a bug.
pdflatex -output-directory="$op_dir" "$op_tex"
pdflatex -output-directory="$op_dir" "$op_tex"
