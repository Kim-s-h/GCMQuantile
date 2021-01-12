#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch
#SBATCH -J N500
#SBATCH -t 1:00:00
#SBATCH --mem=4G
#SBATCH -o ../Routput/N500/results_%a.out

module load R/4.0.0-foss-2019b
filename1="../data/N500/gcm_data_${SLURM_ARRAY_TASK_ID}.rda"
filename2="../results/N500/estimates_gcmQ2_${SLURM_ARRAY_TASK_ID}.rda"
tau1="0.5"
tau2="0.5"
echo $filename1
echo $tau1
echo $tau2
Rscript stan_sapelo_gcm_quantile_random.R  ${SLURM_CPUS_PER_TASK} $filename1 $filename2 $tau1 $tau2

