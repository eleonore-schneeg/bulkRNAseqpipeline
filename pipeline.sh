#!/bin/bash

#added to git

# Mini description: Processing the raw sequence data from the initial QC to genomic mapping and post-alignment QC [raw_fastq  > QC > mapping > indexing > QC]
#########################################################
# QC is done with fastqc, SortMeRna
# Mapping is done with STAR
# Indexing and post-alignment QC is done using samtools, RSeQC and picard
#################################################################

# Script definition steps:


#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -o /rds/general/user/nfancy/home/log/pipeline.^array_index^.out
#PBS -e /rds/general/user/nfancy/home/log/pipeline.^array_index^.err
#PBS -J 1-3


cd $PBS_O_WORKIR

# Load the modules

module load fastqc/0.11.5
module load anaconda3/personal
module load picard/2.6.0
module load java/jdk-8u144
module load R
module load samtools/1.3.1
module load star/2.7.1a
SORTMERNAPATH=/rds/general/user/nfancy/home/packages/sortmerna-2.1b

# Definition of directories
inDir=/rds/general/user/nfancy/ephemeral/maria/HCLK2BBXY/IGFQ000853_Owen_26-7-19_3endRNA-seq #change 1
outDirName=$(basename $inDir | awk -F"_" '{print $1"_out"}')
outDir=$(dirname $inDir)/$outDirName
initialQCdir=$outDir/initialQC
outBAMdir=$outDir/BAM
outReadDistributionDir=$outDir/readDistribution
outGeneBodyCoverageDir=$outDir/geneBodyCoverage
outrRNApercentDir=$outDir/rRNApercent

input=$(head -n $PBS_ARRAY_INDEX $inDir/directories.txt | tail -n 1 ) #change 1
prefix=$(basename $input)
U="_"

mkdir -p $outDir
mkdir -p $initialQCdir
mkdir -p $outBAMdir
mkdir -p $outReadDistributionDir
mkdir -p $outGeneBodyCoverageDir
mkdir -p $outrRNApercentDir


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	InitialQC with fastqc
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo -e "STEP1 started: InitialQC at $(date) "
fastqc -o $initialQCdir -d $TMPDIR $input/*R1_001.fastq.gz
fastqc -o $initialQCdir -d $TMPDIR $input/*R2_001.fastq.gz
echo -e "STEP1 ended: InitialQC at $(date)"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Mapping with STAR and indexing BAM to Bai using samtools
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo -e "STEP2 started:  Mapping reads to hg38 and giving output as BAM at $(date)"

refGenome=/rds/general/user/nfancy/home/STAR_genome/genome.hg38.ucsc.ercc.star.2.7.1a

for i in $input/*R1_001.fastq.gz
do
name=$(basename $i)
name=${name%_R1_001.fastq.gz}

STAR --genomeDir $refGenome --genomeLoad LoadAndKeep --readFilesIn $i --readFilesCommand zcat --outFileNamePrefix $outBAMdir/$prefix$U$name$U --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --runThreadN 8 --limitBAMsortRAM 36000000000

rm -r $outBAMdir/$prefix$U$name$U*STARtmp/

samtools index $outBAMdir/$prefix$U$name$U*.bam

echo -e "STEP4 ended:  Mapping reads to hg38 and giving output as BAM at $(date)"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Percent bases mapped to features using picard
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo -e "STEP3 started: Percent bases mapped to features at $(date)"

bamFile=$(ls $outBAMdir/$prefix$U$name${U}*.bam)

refFlat=/rds/general/user/nfancy/home/STAR_genome/ucsc.hg38.refFlat.txt

picard CollectRnaSeqMetrics TMP_DIR=$TMPDIR I=$bamFile O=$outReadDistributionDir/$prefix$U$name.txt REF_FLAT=$refFlat STRAND=FIRST_READ_TRANSCRIPTION_STRAND CHART_OUTPUT=$outReadDistributionDir/$prefix$U$name.pdf

echo -e "STEP3 ended: Percent bases mapped to features at $(date)"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Genebody coverage using RSeqC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo -e "STEP4 started: Genebody coverage at $(date)"

refBED=/rds/general/user/nfancy/home/STAR_genome/gencode.v29.annotation.bed

geneBody_coverage.py  -r $refBED -i $outBAMdir/$prefix$U$name$U*.bam -o  $outGeneBodyCoverageDir/$prefix$U$name

echo -e "STEP4 ended: Genebody coverage at $(date)"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	rRNA percentage by SortMeRNA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo -e "STEP5 started: rRNA percentage quantification at $(date)"

readFastq=${i%.gz}

gunzip -c $i > $readFastq


$SORTMERNAPATH/sortmerna --ref $SORTMERNAPATH/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNAPATH/index/silva-euk-18s-db:$SORTMERNAPATH/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNAPATH/index/silva-euk-28s-db:$SORTMERNAPATH/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNAPATH/index/rfam-5s-db:$SORTMERNAPATH/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNAPATH/index/rfam-5.8s-db --reads $readFastq --fastx --aligned $outrRNApercentDir/$prefix$U$name.txt --num_alignments 1 --log -m 4096 -v

echo -e "STEP5 ended: rRNA percentage quantification at $(date)"

done
