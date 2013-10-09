#!/bin/sh

# this is a quick demo on how to use miraEST for clean assembly 
#  of ESTs from different strains and subsequent SNP search
#

ln -f -s ../data/fasta_estset1/tvc_mini.fasta tvc_in.sanger.fasta
ln -f -s ../data/fasta_estset1/tvc_mini.fasta.qual tvc_in.sanger.fasta.qual

echo "Running miraSearchESTSNPs, step 1"
miraSearchESTSNPs --project=tvc --job=denovo,normal,sanger,esps1 >log_esps1.txt
echo "Running miraSearchESTSNPs, step 2"
miraSearchESTSNPs --job=denovo,normal,sanger,esps2 >log_esps2.txt
echo "Running miraSearchESTSNPs, step 3"
miraSearchESTSNPs --job=denovo,normal,sanger,esps3 >log_esps3.txt
echo "Done."
echo "If all went well, the results are in 'step3_assembly' directory."

echo "Have a look at the different tag types set in the step3_results files"
echo " how they correspond to the simulated straindata in"
echo " step1_straindata_in.txt"
