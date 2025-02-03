#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=Align-parallel
#SBATCH --partition=long
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk
#SBATCH --output=logs/Align-parallel_%A_%a.log
#SBATCH --error=logs/Align-parallel_%A_%a.error
#SBATCH --array=0-9

# Load modules
ml SAMtools/1.16.1-GCC-11.3.0
export PATH=/data/biol-silvereye/ball6625/software/bwa-mem2-2.2.1_x64-linux:$PATH

# Paths
REF=/data/biol-silvereye/ref_genome/Zlat_2_Tgut_pseudochromosomes.shortChromNames.fasta.gz
RAW_READS_PATH=/data/biol-silvereye/ball6625/norfolk-island/outputs/filtered-reads
OUT_PATH=/data/biol-silvereye/ball6625/norfolk-island/outputs/bams
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# Begin alignment
echo "Aligning $SAMPLE..."

FORWARD=${RAW_READS_PATH}/${SAMPLE}_trimmed_R1.fastq.gz
REVERSE=${RAW_READS_PATH}/${SAMPLE}_trimmed_R2.fastq.gz

# Run bwa-mem2 and directly pipe to SAMtools
bwa-mem2 mem -t 4 $REF $FORWARD $REVERSE | \
samtools view -@ 4 -bS - | \
samtools sort -@ 4 -o ${OUT_PATH}/${SAMPLE}.sorted.bam

# Stats
samtools stats ${OUT_PATH}/${SAMPLE}.sorted.bam > ${OUT_PATH}/${SAMPLE}.stats

# Index bam
samtools index ${OUT_PATH}/${SAMPLE}.sorted.bam

echo "$SAMPLE alignment complete."