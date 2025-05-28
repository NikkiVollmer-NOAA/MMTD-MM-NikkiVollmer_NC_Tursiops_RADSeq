#!/bin/bash
#SBATCH --job-name STRUCT
#SBATCH --mail-user=nicole.vollmer@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/nvollmer/structure/log
#SBATCH --array=[1-100]%25
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=7-00

module load bio/structure/2.3.4

cd ~/structure

# Read line $SLURM_ARRAY_TASK_ID from tasks.txt
K=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks.txt | cut -f1)
RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks.txt | cut -f2)

structure -K $K -m mainparams.txt -o k${K}_run${RUN} 2>&1 | tee k${K}_run${RUN}.log
