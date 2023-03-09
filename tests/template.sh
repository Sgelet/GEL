#!/bin/sh
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J Skeletonization
### -- ask for number of cores (default: 1) --
#BSUB -n 8
### -- specify that the cores must be on the same host --
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot --
#BSUB -R "rusage[mem=1GB]"
### -- Request specific CPU
#BSUB -R "select[model == XeonGold6226R]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot --
#BSUB -M 2GB
### -- set walltime limit: hh:mm --
#BSUB -W 1:00
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
export CC=gcc

# Run
echo FILENAME | ./runtime_test.sh 3
