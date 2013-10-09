#!/bin/sh

mkdir origdata
cd origdata
echo "Fetching sequence of reference genome Escherichia_coli_B_REL606 at GenBank"
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_B_REL606/NC_012967.gbk

echo "Fetching sequencing data (2 files) for strain REL8593A"
wget ftp://ftp.ncbi.nlm.nih.gov/sra/static/SRX012/SRX012992/SRR030257_1.fastq.gz
wget ftp://ftp.ncbi.nlm.nih.gov/sra/static/SRX012/SRX012992/SRR030257_2.fastq.gz

cd ..
