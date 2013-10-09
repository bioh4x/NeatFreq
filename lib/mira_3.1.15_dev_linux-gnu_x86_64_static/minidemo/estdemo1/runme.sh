#!/bin/sh

# this is a quick demo on how to use miraEST for EST assembly and
#  SNP analysis
#
# For the time being, the input files MUST have the prefix 'step1_in' to 
#  be processed correctly by miraeST, I'm sorry.

echo "This demo must be reworked wor 2.9.x and is currently not available, sorry."
exit

echo "Cleaning the provided fasta file"
uncover_at ../data/fasta_estset1/triphysaria_versicolor.masked.fasta -P "-CL:mbc=on:mbcmeg=240" -n -a -r ../data/fasta_estset1/triphysaria_versicolor.fasta step1_in.fasta >ttt
ln -s ../data/fasta_estset1/triphysaria_versicolor.fasta.qual step1_in.fasta.qual
touch step1_straindata_in.txt

echo "Running miraSearchESTSNPs"
miraSearchESTSNPs --fasta | tee run.log
echo "Done."
echo "If all went well, the results are in 'step3_results' directory."
echo "Consult run.log and/or stderr if errors occured."
