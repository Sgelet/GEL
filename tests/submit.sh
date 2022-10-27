#!/bin/sh
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J Skeletonization
### -- ask for number of cores (default: 1) --
#BSUB -n 16
### -- specify that the cores must be on the same host --
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot --
#BSUB -R "rusage[mem=1GB]"
### -- Request specific CPU
##BSUB -R "select[model == XeonGold6226R]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot --
#BSUB -M 3GB
### -- set walltime limit: hh:mm --
#BSUB -W 01:00
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o Output_%J.out
#BSUB -e Error_%J.err

# here follow the commands you want to execute
module load gcc/12.1.0-binutils-2.38 ninja/1.10.2 cmake/3.23.2
rm -rf ../cmake-build-default
rm -f out_raw/*
export CC=gcc
unzip -o 3dmeshes.zip -d /work3/etoga/

# Variations to run
declare -a variations=(
  ""
)

variant=0

mkdir -p out_raw
mkdir -p skeletons

for i in "${variations[@]}"
do
  # Build with compile time option
  cmake -G Ninja -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_CXX_FLAGS="$i -O2 -DCORE_TEST=$LSB_DJOB_NUMPROC -DCORE_TEST_SEC=2" -S .. -B ../cmake-build-default -DCMAKE_BUILD_TYPE=Release
  cmake --build ../cmake-build-default

  # Make directories
  mkdir -p "skeletons/var$((++variant))"

  # Run
  ls -1 /work3/etoga/3DMeshes/arm* | ./runtime_test.sh 3 "var$variant"
done

# Cleanup
rm -f /work3/etoga/3DMeshes/*