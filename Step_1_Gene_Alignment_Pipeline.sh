#!/bin/bash
#Author: Anurag Kanase
#Github: dezember
#Version: 1
#Date: August 12, 2020
#Project: Precision Drug Prediction
echo "Welcome to Precision Drug Prediction"
echo "Running Step 1 of automated Pipeline"
echo "Download dependencies requirements"
#sudo apt install subread parallel samtools
conda install -c bioconda bowtie2 sra-tools parallel-fastq-dump
sudo apt get samtools

echo "Creating Directories"
mkdir /home/anurag/Seq
mkdir /home/anurag/Seq/Analysis
mkdir /home/anurag/Seq/Healthy
mkdir /home/anurag/Seq/Infected
mkdir /home/anurag/Seq/SAM
mkdir /home/anurag/Seq/BAM
mkdir /home/anurag/Seq/ref
mkdir /home/anurag/Seq/sorted_BAM
mkdir /home/anurag/Seq/GFF3
mkdir /home/anurag/Seq/Counts 

echo "Creating variable file names. H = Healthy, D = Infected/Diseased"

cd Seq
H1="SRR11804725"
H2="SRR11804726"
H3="SRR11804727"
H4="SRR11804728"

D1="SRR11804721"
D2="SRR11804722"
D3="SRR11804723"
D4="SRR11804724"

echo "Starting SRA File Downloading. This process will vary on download time"
echo "Downloading Control Sequences"
prefetch ${H1}
prefetch ${H2}
prefetch ${H3}
prefetch ${H4}
echo "Downloading Diseased Sequences"
prefetch ${D1}
prefetch ${D2}
prefetch ${D3}
prefetch ${D4}

echo "SRR Download Finished."

echo "Converting SRA to FASTQ files."
#Pair END
parallel-fastq-dump --sra-id ${H1} --threads 16 --outdir Healthy/ --split-files --gzip
parallel-fastq-dump --sra-id ${H2} --threads 16 --outdir Healthy/ --split-files --gzip
parallel-fastq-dump --sra-id ${H3} --threads 16 --outdir Healthy/ --split-files --gzip
parallel-fastq-dump --sra-id ${H4} --threads 16 --outdir Healthy/ --split-files --gzip

echo "Healthy FASTQ dumping over"
parallel-fastq-dump --sra-id ${D1} --threads 16 --outdir Infected/ --split-files --gzip
parallel-fastq-dump --sra-id ${D2} --threads 16 --outdir Infected/ --split-files --gzip
parallel-fastq-dump --sra-id ${D3} --threads 16 --outdir Infected/ --split-files --gzip
parallel-fastq-dump --sra-id ${D4} --threads 16 --outdir Infected/ --split-files --gzip


echo "Downloading Reference Sequence"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
tar -xzvf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index* ref/

echo "Aligning the sequenced data to reference genome"


#for-pair-end
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Healthy/${H1}_1.fastq.gz -2 Healthy/${H1}_2.fastq.gz -S SAM/${H1}.sam
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Healthy/${H2}_1.fastq.gz -2 Healthy/${H2}_2.fastq.gz -S SAM/${H2}.sam
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Healthy/${H3}_1.fastq.gz -2 Healthy/${H3}_2.fastq.gz -S SAM/${H3}.sam
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Healthy/${H4}_1.fastq.gz -2 Healthy/${H4}_2.fastq.gz -S SAM/${H4}.sam
echo "Healthy bt2 sequenced finished"

bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Infected/${D1}_1.fastq.gz -2 Infected/${D1}_2.fastq.gz -S SAM/${D1}.sam
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Infected/${D2}_1.fastq.gz -2 Infected/${D2}_2.fastq.gz -S SAM/${D2}.sam
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Infected/${D3}_1.fastq.gz -2 Infected/${D3}_2.fastq.gz -S SAM/${D3}.sam
bowtie2 -p 16 -x ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 Infected/${D4}_1.fastq.gz -2 Infected/${D4}_2.fastq.gz -S SAM/${D4}.sam

echo "Converting SAM to BAM"
samtools view -b -o BAM/${H1}.bam SAM/${H1}.sam
samtools view -b -o BAM/${H2}.bam SAM/${H2}.sam
samtools view -b -o BAM/${H3}.bam SAM/${H3}.sam
samtools view -b -o BAM/${H4}.bam SAM/${H4}.sam

samtools view -b -o BAM/${D1}.bam SAM/${D1}.sam
samtools view -b -o BAM/${D2}.bam SAM/${D2}.sam
samtools view -b -o BAM/${D3}.bam SAM/${D3}.sam
samtools view -b -o BAM/${D4}.bam SAM/${D4}.sam

echo "Sorting BAM files"
samtools sort -@ 16 -m 3G -o sorted_BAM/${H1}.sort.bam BAM/${H1}.bam
samtools sort -@ 16 -m 3G -o sorted_BAM/${H2}.sort.bam BAM/${H2}.bam
samtools sort -@ 16 -m 3G -o sorted_BAM/${H3}.sort.bam BAM/${H3}.bam
samtools sort -@ 16 -m 3G -o sorted_BAM/${H4}.sort.bam BAM/${H4}.bam

samtools sort -@ 16 -m 3G -o sorted_BAM/${D1}.sort.bam BAM/${D1}.bam
samtools sort -@ 16 -m 3G -o sorted_BAM/${D2}.sort.bam BAM/${D2}.bam
samtools sort -@ 16 -m 3G -o sorted_BAM/${D3}.sort.bam BAM/${D3}.bam
samtools sort -@ 16 -m 3G -o sorted_BAM/${D4}.sort.bam BAM/${D4}.bam

echo "Downloading annotation file"
#extracting 
wget ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz
gunzip Homo_sapiens.GRCh38.100.gff3.gz
mv Homo_sapiens.GRCh38.100.gff3 /home/anurag/Seq/GFF3

#Starting Feature Counting
echo "Merging Differential Tables"
ANNOT_GFF="/home/anurag/Seq/GFF3/Homo_sapiens.GRCh38.100.gff3"
INDIR="/home/anurag/Seq/BAM"
OUTDIR="/home/anurag/Seq/Counts"

parallel -j 4 "featureCounts -T 4 -s 2 -p -t gene -g ID -a ${ANNOT_GFF} -o ${OUTDIR}/{/.}.gene.txt {}" ::: ${INDIR}/*.bam

echo "Counting Features for differential expression"
python /home/anurag/features.py

#Converting ENSEMBL names to GENE names
pip install pandas
pip install pyensembl
pyensembl install --release 100 --species homo_sapiens