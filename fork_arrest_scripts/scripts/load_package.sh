#!/bin/bash

# script that loads packages on the hpc
# usage: source load_package.sh -name1 -name2 -name3
# etc.

# To use the programs in the repo, you need the pre-requisite software packages.
# There are two ways to get the packages.
# (1) You can make a singularity image using the singularity.def file we have provided
#     with the repo. If you do this, please edit the singularity_img variable below
#     to point to this image. If you follow this method, please read the comments in the
#     singularity.def file first.
# (2) You've obtained the names of packages from the singularity.def file but have installed
#     all the packages directly on your system. In this case,
#     please delete all the lines from this file but do not delete the file itself.
#     This will ensure lines like `source load_package.sh -python` in our scripts
#     do not do anything but do not throw an error either. NOTE: for historical reasons,
#     a small number of our packages below like paf2chain have to be called as $paf2chain
#     instead of just paf2chain. So if you have used method (2), please ensure that these
#     packages can be called using this $ notation i.e. that their invocation is stored
#     in a shell variable.

# You will see ":" in the code below. This is a no-op command in bash.
# Some of the packages below are not available in the singularity image we provide, due to different reasons
# such as license issues, future-proofing this repo so they are not needed right now etc. You are free to
# install them on your system and use them. If you do so, please fix these invocations.
# Also see the comments in the singularity.def file for more details.

# supply the path below
singularity_img=/path/to/singularity/image.img; # this may also end in a .sif

export paftools_js="singularity exec $singularity_img paftools.js "
export paf2chain="singularity exec $singularity_img paf2chain "

python() {
  singularity exec "$singularity_img" python3.10 "$@";
}
samtools() {
  singularity exec "$singularity_img" samtools "$@";
}
seqtk() {
  singularity exec "$singularity_img" seqtk "$@";
}
minimap2() {
  singularity exec "$singularity_img" minimap2 "$@";
}
bedtools() {
  singularity exec "$singularity_img" bedtools "$@";
}
pycoQC() {
  singularity exec "$singularity_img" pycoQC "$@";
}
R() {
  singularity exec "$singularity_img" R --no-save "$@";
}
Rscript() {
  singularity exec "$singularity_img" Rscript "$@";
}
pdflatex() {
  singularity exec "$singularity_img" pdflatex "$@";
}
pdflatex_w_tikz() {
  singularity exec "$singularity_img" pdflatex "$@";
}
convert() {
  singularity exec "$singularity_img" convert "$@";
}
mlr() {
  singularity exec "$singularity_img" mlr "$@";
}
jq() {
  singularity exec "$singularity_img" jq "$@";
}
modkit() {
  singularity exec "$singularity_img" modkit "$@";
}
seqkit() {
  singularity exec "$singularity_img" seqkit "$@";
}
bedops() {
  singularity exec "$singularity_img" bedops "$@";
}
wig2bed() {
  singularity exec "$singularity_img" wig2bed "$@";
}
git() {
  singularity exec "$singularity_img" git "$@";
}

for var in "$@"
do
    if [ "$var" == "-python" ]
    then
      export -f python;
    elif [ "$var" == "-samtools" ]
    then
        export -f samtools;
    elif [ "$var" == "-seqtk" ]
    then
        export -f seqtk;
    elif [ "$var" == "-guppy_5_07" ]
    then
        :
    elif [ "$var" == "-guppy_6_57" ]
    then
        :
    elif [ "$var" == "-git" ]
    then
        export -f git;
    elif [ "$var" == "-minimap2pt24" ] || [ "$var" == "-paftools_js" ]
    then
        export -f minimap2;
        # we've already defined paftools_js earlier in this script
    elif [ "$var" == "-dnascent" ]
    then
        :
    elif [ "$var" == "-dnascent_v4" ]
    then
        :
    elif [ "$var" == "-bedtools" ]
    then
        export -f bedtools;
    elif [ "$var" == "-pycoQC" ]
    then
        export -f pycoQC;
    elif [ "$var" == "-R" ]
    then
        export -f R;
        export -f Rscript;
    elif [ "$var" == "-latex" ]
    then
      export -f pdflatex;
    elif [ "$var" == "-pdflatex_w_tikz" ]
    then
      export -f pdflatex_w_tikz;
    elif [ "$var" == "-imagemagick" ]
    then
      export -f convert;
    elif [ "$var" == "-miller" ]
    then
      export -f mlr;
    elif [ "$var" == "-jq" ]
    then
      export -f jq;
    elif [ "$var" == "-bigWigToBedGraph" ]
    then
      :
    elif [ "$var" == "-dorado_0_3_4" ]
    then
      :
    elif [ "$var" == "-dorado_0_6_2" ]
    then
      :
    elif [ "$var" == "-dorado_0_7_2" ]
    then
      :
    elif [ "$var" == "-remora_3_1_0" ]
    then
      :
    elif [ "$var" == "-remora_3_2_0" ]
    then
      :
    elif [ "$var" == "-bedops" ]
    then
      export -f bedops;
      export -f wig2bed;
    elif [ "$var" == "-deeptools" ]
    then
      :
    elif [ "$var" == "-modkit" ]
    then
      export -f modkit;
    elif [ "$var" == "-pod5" ]
    then
      :
    elif [ "$var" == "-liftOver" ]
    then
      :
      # we've not provided liftOver in the def file. Please see comments in that file for details.
      # export liftOver="singularity exec something.img liftOver"
    elif [ "$var" == "-seqkit" ]
    then
      export -f seqkit;
    elif [ "$var" == "-paf2chain" ]
    then
      # we've already defined paf2chain earlier in this script
      :
    fi
done
