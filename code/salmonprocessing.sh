#!/bin/bash

conda install fastqc 

fastqc *.fastq

# Trim

conda install fastp

for f in *.fastq
do
    fastp -i $f -o trimmed_$f -h ${f%.fastq}_fastp.html -j ${f%.fastq}_fastp.json
done

# fastqc again

fastqc trimmed_*.fastq

# Worse results, using original data

# Download reference transcriptome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz

gunzip GCF_000146045.2_R64_rna.fna.gz

# Convert reads to gene counts

conda install salmon

salmon index -t GCF_000146045.2_R64_rna.fna -i yeastindex

# Run for all 9

for f in *.fastq
do
  salmon quant -i yeastindex -l A -r $f -p 4 -o ${f%.fastq}_quant
done

# Now have quant.sf files need for R





