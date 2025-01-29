#!/bin/bash
#SBATCH --time=3:00:00            
#SBATCH --job-name=fastp     
#SBATCH --partition=short
#SBATCH --array=1-10
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8        
#SBATCH --output=filter_%A.log      
#SBATCH --error=filter_%A.error
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk

# Define input and output directories
INPUT_DIR=/data/biol-silvereye/norfolk_wgs/arbor
OUTPUT_DIR=/data/biol-silvereye/ball6625/norfolk-island/filter/filtered_reads
mkdir -p $OUTPUT_DIR

# Extract input files for this task
R1=$(ls $INPUT_DIR/*_R1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
R2=$(ls $INPUT_DIR/*_R2.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
OUT_R1=${OUTPUT_DIR}/$(basename "$R1" | sed 's/_R1/_trimmed_R1/')
OUT_R2=${OUTPUT_DIR}/$(basename "$R2" | sed 's/_R2/_trimmed_R2/')

# Run fastp, borrowing settings from Norfolk Island paper
/data/biol-silvereye/ball6625/software/fastp -i "$R1" -I "$R2" -o "$OUT_R1" -O "$OUT_R2" -w 8 \
--trim_front1 10 --trim_front2 10 --detect_adapter_for_pe --dedup
