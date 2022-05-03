#!/usr/bin/env bash

#RNA sequencing analyis pipeline using GRCh38 genome as a reference. Utilizes FASTQC, MultiQC as quality control to
#ensure seequences are up to standards. Utilizes both samtools, hisat2, and featureCounts to run alignment and analysis.
#Optional visualization commands using bedtools and BedGraphToBigWig to view alignments.

#__author__ = "T-Chinsky"
#__copyright__ = "Copyright 2022"
#__credits__ = ["T-Chinsky"]
#__license__ = "GNU V3"
#__version__ = "1.0.0"
#__maintainer__ = "T-Chinsky"
#__status__ = "Completed"

set -e

#prints out last command when errors occur
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "\"${last_command}\" command filed with code $?."' EXIT

#Getting reference genome
cd refs/

wget -O GRCh38.p13_genomic.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
wget -O GRCh38.p13_genomic.gtf.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz"
wget -O GRCh38.p13_rna.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz"

gunzip *.gz

cd ..

# building the reference genome
echo "Building reference genome..."

IDX=refs/GRCh38.p13_genomic.fna
printf "###############################################                                                  (50%%)\n"
hisat2-build -p 32 $IDX $IDX
printf "################################################################################################ (100%%)\n"
printf '\n'

#Index the reference genome with samtools
echo "Indexing reference genome..."

printf "###############################################                                                  (50%%)\n"
samtools faidx $IDX
printf "################################################################################################ (100%%)\n"
printf '\n'

#Running QC on all files
echo "Running inital QC on FASTQ files..."
mkdir female_raw_qc
mkdir male_raw_qc
printf "################                                                                                 (20%%)\n"
cat female_ids | parallel "fastqc reads/{}_R1.fastq -o female_raw_qc/"
printf "#####################################                                                            (40%%)\n"
cat female_ids | parallel "fastqc reads/{}_R2.fastq -o female_raw_qc/"
printf "#########################################################                                        (60%%)\n"
cat male_ids | parallel "fastqc reads/{}_R1.fastq -o male_raw_qc/"
printf "#############################################################################                    (80%%)\n"
cat male_ids | parallel "fastqc reads/{}_R2.fastq -o male_raw_qc/"
printf "################################################################################################ (100%%)\n"
printf '\n'
echo "Generating QC reports..."
mkdir female_raw_multiqc
mkdir male_raw_multiqc
printf "######                                                                                           (10%%)\n"
multiqc female_raw_qc/ -o female_raw_multiqc/
printf "###############################################                                                  (50%%)\n"
multiqc male_raw_qc/ -o male_raw_multiqc/
printf "################################################################################################ (100%%)\n"
printf '\n'

#Store index name in variable
IDX=refs/GRCh38.p13_genomic.fna

#Create the bam folder for output
echo "Creating output folders for bam files..."

mkdir -p female_bam
mkdir -p male_bam

#Align the FASTQ files to the reference genome
echo "Aligning female FASTQ files to the reference genome..."

cat female_ids | parallel "hisat2 -x $IDX -1 reads/{}_R1.fastq -2 reads/{}_R2.fastq | samtools sort > female_bam/{}.bam"

printf "###############################################                                                  (50%%)\n"
echo "Aligning male FASTQ files to the reference genome..."

cat male_ids | parallel "hisat2 -x $IDX -1 reads/{}_R1.fastq -2 reads/{}_R2.fastq | samtools sort > male_bam/{}.bam"
printf "################################################################################################ (100%%)\n"
printf '\n'

#Index each BAM file
echo "Indexing the bam files..."

cat female_ids | parallel "samtools index female_bam/{}.bam"
printf "###############################################                                                  (50%%)\n"
cat male_ids | parallel "samtools index male_bam/{}.bam"
printf "################################################################################################ (100%%)\n"
printf '\n'

#Run QC on Aligned files
echo "Running alignment QC..."
mkdir male_align_qc
mkdir female_align_qc
printf "##############################                                                                   (33%%)\n"
cat female_ids | parallel "fastqc female_bam/{}.bam -o female_align_QC/"
printf "#############################################################                                    (66%%)\n"
cat male_ids | parallel "fastqc male_bam/{}.bam -o male_align_QC/"
printf "################################################################################################ (100%%)\n"
printf '\n'
echo "Generating reports..."
mkdir female_align_multiqc
mkdir male_align_multiqc

multiqc female_align_qc -o female_align_multiqc/
printf "###############################################                                                  (50%%)\n"
multiqc male_align_qc -o male_align_multiqc/
printf "################################################################################################ (100%%)\n"
printf '\n'


##Optional bedcoverage and bigwig steps to visualize the alignments
#Convert bam files to bedgraph coverage files
#cat female_ids | parallel "bedtools genomecov -ibam female_bam/{}.bam -split -bg > female_bam/{}.bg"
#cat male_ids | parallel "bedtools genomecov -ibam male_bam/{}.bam -split -bg > male_bam/{}.bg"

#Convert each bedgraph coverage into bigwig coverage
#cat female_ids | parallel "bedGraphToBigWig female_bam/{}.bg  ${IDX}.fai female_bam/{}.bw"
#cat male_ids | parallel "bedGraphToBigWig male_bam/{}.bg  ${IDX}.fai male_bam/{}.bw"

#Get feature counts summarizing reads
echo "Getting feature counts for each file..."

featureCounts -p -T 32 -a refs/GRCh38.p13_genomic.gtf -o female_counts.txt female_bam/CR???.bam female_bam/PV???.bam female_bam/PV???_?.bam
printf "###############################################                                                  (50%%)\n"
featureCounts -p -T 32 -a refs/GRCh38.p13_genomic.gtf -o male_counts.txt male_bam/CR???.bam male_bam/PV???.bam male_bam/PV???_?.bam
printf "################################################################################################ (100%%)\n"
printf '\n'

echo "Completed Alignment analysis..."
