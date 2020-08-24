# CodeforCovid19_hackathon
We have developed a very simple three-step pipeline:
1. to find the NGS SRA files from NCBI portal of healthy vs diseased patients, align it to the Human Genome GRCh38.pt13 reference sequence using Bowtie-2 alignment, and find the featureCounts of the genes in individual files. 
2. in the second step we use the features_1.csv generated in #1. to find out the differential gene expressions using R. 
3. The differential gene expression csv file generated in #2. will be used for automatically find the gene names from their existing Ensembl names, and directly find the drug target associated with them. 

The Drug discovery analysis can be performed as follows:

### 1. Gene alignment
Run the line below in shell:
<pre><code>sh .\Step_1_Gene_Alignment_Pipeline.sh</pre></code>
Note: you make need to edit the number of cores used depending on the processing power of your computer. We have currently use -16 cores. This pipeline is the longest time consumer. If you are going to use whole human genome for analysis, the average time may vary between 8 hours to 36 hours for 24 cores to 4 cores. 

### 2. Differential Gene Expressions
Run the file in R, using the GUI of RStudio or RScript. The prerequisite of Bioconductor and DESeq2 is already included in the file. 

### 3. Drug Discovery
Run the file using Anaconda or Python environment. 

<pre><code> python .\Step_3_Finding_Drugs.py </pre></code>

For step #1, it's important to have Linux environment. GPU is not necessary. Depending on the type of analysis, 1 TB free space and minimum 30 GB RAM is recommended. 

Step #2 and Step #3 can be performed on any OS. 
