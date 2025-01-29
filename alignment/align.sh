#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --job-name=Align
#SBATCH --partition=long
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk
#SBATCH --output=Align_%A.log
#SBATCH --error=Align_%A.error

# Load modules
ml SAMtools/1.16.1-GCC-11.3.0

# Define paths
export PATH=/data/biol-silvereye/ball6625/software/bwa-mem2-2.2.1_x64-linux:$PATH
REF=/data/biol-silvereye/ref_genome/GCA_001281735.1_ASM128173v1_genomic.fna.gz
RAW_READS_PATH=/data/biol-silvereye/ball6625/norfolk-island/filter/filtered-reads
mapfile -t SAMPLES < samples.txt

# Index reference
if [ ! -f "$REF.bwt.2bit.64" ]; then
    echo "Indexing reference genome..."
    bwa-mem2 index $REF
fi

#/data/biol-silvereye/ball6625/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index $REF

# Proceed to alignment

for SAMPLE in "${SAMPLES[@]}"
do
    echo "Aligning $SAMPLE..."

    # point to forward and reverse files
    FORWARD=${RAW_READS_PATH}/${SAMPLE}_trimmed_R1.fastq.gz
    REVERSE=${RAW_READS_PATH}/${SAMPLE}_trimmed_R2.fastq.gz

    # Run bwa-mem2 for PE reads
    #/data/biol-silvereye/ball6625/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem \
    #    -t 8 \
    #    $REF \
    #    $FORWARD \
    #    $REVERSE \
    #    > ${SAMPLE}.sam
    
    bwa-mem2 mem -t 8 $REF $FORWARD $REVERSE > ${SAMPLE}.sam

    # check stats
    samtools stats ${SAMPLE}.sam > ${SAMPLE}.stats

    # convert sam to bam
    samtools view -bS ${SAMPLE}.sam > ${SAMPLE}.bam

    # sort bam
    samtools sort ${SAMPLE}.bam -o ${SAMPLE}.sorted.bam

    # index sorted bam
    samtools index ${SAMPLE}.sorted.bam

    echo "$SAMPLE alignment complete."

done

echo "All alignments complete."

