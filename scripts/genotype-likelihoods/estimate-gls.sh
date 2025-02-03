#!/bin/bash
#SBATCH --array=1-33
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --job-name=EstimateGLs
#SBATCH --partition=long
#SBATCH --output=logs/EstimateGLs_%A_%a.log
#SBATCH --error=logs/EstimateGLs_%A_%a.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk

#########################################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION WITH ANGSD v.0.921
# Will output the following files per pseudo-chromosome:
# 1. VCF with genotype probability (GP) and genotype likelihood (GL) fields
# 2. MAFs
# 3. Genotype likelihoods (GLs) in beagle format

# A. Sendell-Price, April 2022
# Adapted by Abby Williams, January 2025
#########################################################################################################

#Load angsd module
#ml angsd/0.925-foss-2018b

# Set path to reference assembly, ancestral sequence, bam files and output
# Note: bam files need to be indexed (using samtools index) 
REF=/data/biol-silvereye/ref_genome/Zlat_2_Tgut_pseudochromosomes.shortChromNames.fasta.gz
ANC=/data/biol-silvereye/ball6625/norfolk-island/outputs/outgroup-consensus/full_15-179.fa.gz
OUTPUTS=/data/biol-silvereye/ball6625/norfolk-island/outputs/genotype-likelihoods
BAMS=/data/biol-silvereye/ball6625/norfolk-island/scripts/genotype-likelihoods/bams.list

# Use slurm array task ID to get long chromosome name
CHROM=$(cat chrom.list | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#Make directory for focal chromosome and move into it
OUT_PATH=$OUTPUTS/$CHROM
mkdir -p $OUT_PATH
cd $OUT_PATH

# Estimate genotype likelihoods and output SNPs using ANGSD
angsd -b ${BAMS} -ref $REF -anc $ANC  \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 25 -minQ 25 \
-GL 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGlf 2 -doVcf 1 \
-minInd 5 -minMaf 0.05 -skipTriallelic \
-r $CHROM -out ${CHROM}_polarised_minMapQ25_minQ25_minMAF0.1_minInd5

# Explanation of above settings:
# ==============================
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -GL = calculate genotype likelihoods (1: using SAMtools model )
# -doMajorMinor = infer major and minor alleles (5: Use ancestral allele as major, requires -anc)
# -doMAF = estimate major and minor allele frequencies (1: Known major, and Known minor)
# -doPost = calculate posterior prob (1: Using frequency as prior)
# -doGlf = output the log genotype likelihoods to a file (2: beagle genotype likelihood format)
# -doVcf = output a VCF file (1: yes)
# -minInd = only output sites with information for at least [int] samples (half is sensible)
# -minMaf = Remove sites with MAF below [float]
# -SNP_pval = Remove sites with a pvalue larger [float]
# -r = region to output (chr name e.g. chr1)
# -out = prefix to use for output files
# ==============================

#Dsuite (which we will use for calculating D-stats) can only take VCFs as input but has trouble handling
#the VCF file directly outputted by angsd. To overcome this issue we will use genotype imputation software
#beagle to convert the angsd genotype likelihood format file to a dsuite friendly VCF file.
#Note: the outputted vcf will contain imputed genotypes (GT field), dosages - a linear transformation of the posterior
#genotype probabilities (DS field), and genotype probabilities (GP field). The imputed genotypes will be disregarded
#as we will run Dsuite using genotype likelihood mode.

java -Xmx15000m -jar /data/zool-zost/BIN/beagle.r1399.jar \
gl=${CHROM}_polarised_minMapQ25_minQ25_minMAF0.1_minInd5.vcf.gz \
out=${CHROM}_polarised_minMapQ25_minQ25_minMAF0.1_minInd5.BeagleVersion 
