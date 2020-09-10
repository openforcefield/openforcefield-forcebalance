#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J benchmar
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  %J.o
#BSUB -e  %J.e
#
#BSUB -q cpuqueue

# Set openeye license
export OE_LICENSE="~/oe_license.txt"

# Enable conda
. ~/miniconda3/etc/profile.d/conda.sh

# Use the right conda environment
conda activate propertyestimator

# Launch my program.
module load cuda/9.2

python run.py &> console_output.log
