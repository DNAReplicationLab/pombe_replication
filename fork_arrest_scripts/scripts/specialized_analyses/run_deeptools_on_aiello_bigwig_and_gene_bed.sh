#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J cracDataPerGeneDeepTools
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Do deeptools on CRAC data in 10.1016/j.molcel.2022.06.021 in bigwig format using locations of genes.

# usage
# -----
# sbatch run_deeptools_on_aiello_bigwig_and_gene_bed.sh $prefix_bigwig $ref_fai $gene_bed $tRNA_rDNA_bed $strand \
#   $op_dir [$chrM_optional] [$bin_size] [$average_style]
# NOTE: can use bash instead of sbatch but might take a few minutes.
# NOTE: [] means optional arguments. To specify an optional argument, you must specify all optional arguments to
#       the left of it as well, and you must remove the square brackets for the ones you want to specify.
# NOTE: \ means that the command continues on the next line.
# prefix_bigwig: If files are at /to/file/GSM_R1_minus.bw and ...plus.bw, then prefix is /to/file/GSM_R1.
#                Aiello et al have a plus and a minus bigwig file per experimental condition.
#                So, download those, put them in a folder, and use the prefix corresponding to plus and minus of the
#                same condition here. Also see the strand argument below.
# ref_fai: Fasta index file for reference genome
# gene_bed: Bed file with coding sequence or transcript annotations, must have at least six columns
# tRNA_rDNA_bed: Bed file with tRNA and rDNA annotations, must have at least six columns.
#                Used to remove data in tRNA and rDNA regions on both strands.
# strand: Can be "plus" or "minus" or "". Depending on this input, only that strand of bigwig file
#         i.e. ${prefix_bigwig}_${strand}.bw will be chosen and the gene bed file will be filtered to
#         retain only this strand. But both strands of the exclusion file will be used to perform the exclusion.
#         If set to "", we will not filter the gene bed file by strand, and we will load the bigwig ${prefix_bigwig}.bw
#         i.e. without any strand suffix. This is useful if you want to analyze both strands together.
#         Aiello et al. do not have a "both strand" file, so this option is if you want to do analysis
#         with some other dataset or if you have combined their strand-specific data.
# op_dir: Output directory where the output files will be written.
# chrM_optional: Exclude this mitochondrial chromosome and "chrMito" from the analysis.
#                Default is to exclude "chrM" and "chrMito". We have to include chrMito as that is what Aiello
#                et al call the chromosome. If you are doing some analysis including mitochondrial genes etc.
#                then you may have to manually edit the file suitably in this place (and maybe others).
# bin_size: Bin size for deeptools. Default is 25.
# average_style: Style of averaging within bins. Default is mean. Refer to the parameter --averageTypeBins
#                in the deepTools computeMatrix documentation for more details.
#                (https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html)

# outputs
# -------
# Many outputs from deeptools with names like $op_dir/${base name of prefix_bigwig}_${strand}.mat.gz and
# $op_dir/${base name of prefix_bigwig}_${strand}.png.
# The png files are what we are looking for to put in the paper.

# stop execution if any command fails
set -e

# assign arguments to variables
prefix=${1:-}
ref_fai=${2:-}
gene_bed=${3:-}
tRNA_rDNA_bed=${4:-}
strand=${5:-}
op_dir=${6:-}
chrM=${7:-"chrM"}
bin_size=${8:-25}
average_style=${9:-mean}

# convert paths to absolute paths
prefix=$(realpath "$prefix")
ref_fai=$(realpath "$ref_fai")
gene_bed=$(realpath "$gene_bed")
tRNA_rDNA_bed=$(realpath "$tRNA_rDNA_bed")

# go to the main directory for scripts
cd ..;

# load packages
source load_package.sh -python -deeptools -bigWigToBedGraph

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# check that the correct number of arguments were provided
if [ "$#" -lt 6 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch convert_aiello_bigwig_to_per_gene_bed.sh \$prefix_bigwig \$ref_fai \$gene_bed \$tRNA_rDNA_bed \$strand \$op_dir [\$chrM_optional] [\$bin_size] [\$average_style]"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "Exiting."
    exit 1;
fi

# set up some file names and make the output directory
mkdir -p "$op_dir"
bw_file="$prefix".bw
output_file_prefix="$op_dir"/"$(basename "$prefix")"

# check that strand is either plus or minus, and set the strand sign accordingly
if [ ! "$strand" == "plus" ] && [ ! "$strand" == "minus" ] && [ ! "$strand" == "" ]; then
    >&2 echo "ERROR: strand must be either plus or minus or empty i.e. \"\""
    >&2 echo "Exiting."
    exit 1;
elif [ "$strand" == "plus" ]; then
    strand_sign="+"
    bw_file="$prefix"_"${strand}".bw
    output_file_prefix="$op_dir"/"$(basename "$prefix")"_"${strand}"
elif [ "$strand" == "minus" ]; then
    strand_sign="-"
    bw_file="$prefix"_"${strand}".bw
    output_file_prefix="$op_dir"/"$(basename "$prefix")"_"${strand}"
else
    strand_sign=""
    bw_file="$prefix".bw
    output_file_prefix="$op_dir"/"$(basename "$prefix")"
fi

# check that input files exist
if [ ! -f "$bw_file" ]; then
    >&2 echo "ERROR: $bw_file does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$ref_fai" ]; then
    >&2 echo "ERROR: $ref_fai does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$gene_bed" ]; then
    >&2 echo "ERROR: $gene_bed does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$tRNA_rDNA_bed" ]; then
    >&2 echo "ERROR: $tRNA_rDNA_bed does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

# ensure that the input bed files are valid
if [ ! "$(< "$gene_bed" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: bed file $gene_bed is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$gene_bed" python validate_bed_against_fai.py "$ref_fai" )" == "valid"  ]; then
  >&2 echo "Error: bed file $gene_bed does not have valid coordinates."
  exit 1;
fi

if [ ! "$(< "$tRNA_rDNA_bed" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: bed file $tRNA_rDNA_bed is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$tRNA_rDNA_bed" python validate_bed_against_fai.py "$ref_fai" )" == "valid"  ]; then
  >&2 echo "Error: bed file $tRNA_rDNA_bed does not have valid coordinates."
  exit 1;
fi

# convert input bigwig to bed for the sole purpose of checking if it is valid
tmp_file_bg=$(mktemp "$tmpDir"/tmp_file.XXXXXXXXXX)
tmp_file_bed=$(mktemp "$tmpDir"/tmp_file.XXXXXXXXXX)
bigWigToBedGraph "$bw_file" "$tmp_file_bg"
< "$tmp_file_bg" awk -v OFS="\t" '{print $1,$2,$3,"blank",$4,"+"}' | sort -k 1,1 -k2,2n > "$tmp_file_bed"

if [ ! "$(< "$tmp_file_bed" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: big wig file $bw_file is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$tmp_file_bed" grep -v chrMito | python validate_bed_against_fai.py "$ref_fai" )" == "valid"  ]; then
  >&2 echo "Error: big wig file $bw_file does not have valid coordinates."
  exit 1;
fi

# filter bed file to retain the appropriate strand and only the first six columns and exclude mitochondrial chromosome
tmp_gene_bed=$(mktemp "$tmpDir"/tmp_file.XXXXXXXXXX)
< "$gene_bed" grep -E -v '^browser|^track|^#|^'"$chrM"'|^chrMito' |\
 {
   if [ "$strand_sign" == "" ]; then
     cat
   else
    awk -v OFS="\t" -v strand_sign="$strand_sign" '$6 == strand_sign {print $1,$2,$3,$4,$5,$6}'
  fi
 } > "$tmp_gene_bed"

# perform deeptools analysis
# NOTE: we are sorting by sum as that is how we represent quartiles of genes in the paper
# NOTE: we are plotting median so that we can compare easily with the Aiello et al paper who also
#   plot medians. They do not plot the R loop CRAC data which is why we have to make these plots.
edge_region_bp=0 # we treat regions +- this much near TES, TSS differently
                 # as shown by the options below
computeMatrix scale-regions -S "$bw_file" -R "$tmp_gene_bed" \
                              --beforeRegionStartLength "$edge_region_bp" \
                              --regionBodyLength 1200 \
                              --afterRegionStartLength "$edge_region_bp" \
                              --unscaled5prime "$edge_region_bp" \
                              --unscaled3prime "$edge_region_bp" \
			                        --blackListFileName "$tRNA_rDNA_bed" \
                              --binSize "${bin_size:-25}" \
                              --sortUsing sum \
                              --averageTypeBins "${average_style:-mean}" \
                              --sortRegions descend \
                              --missingDataAsZero -o "$output_file_prefix".mat.gz \
                              --outFileSortedRegions "$output_file_prefix".bed >\
                                "$output_file_prefix".compute_matrix.log 2>&1

for stat in mean median sum; do
  plotHeatmap -m "$output_file_prefix".mat.gz -out "$output_file_prefix"."$stat".png \
    --heatmapHeight 20 --heatmapWidth 10 --colorMap Blues --legendLocation none --dpi 600 \
    --averageTypeSummaryPlot "$stat"
done

# remove temporary directory
rm -rf "$tmpDir"