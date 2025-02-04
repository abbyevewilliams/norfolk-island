#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=short
#SBATCH --time=30:00
#SBATCH --job-name=plot-r

module purge
module load R/4.4.0-gfbf-2023a

Rscript --no-restore --no-save r_script.R