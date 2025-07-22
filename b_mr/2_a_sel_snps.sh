#!/bin/bash
#SBATCH --job-name=clump_proxies
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=3GB
#SBATCH --array=1-2

module add languages/r/4.3.3
R CMD BATCH b_mr/2_a_sel_snps.R b_mr/2_a_sel_snps_$SLURM_ARRAY_TASK_ID.Rout