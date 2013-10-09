#!/bin/sh

# this is a quick demo on how to use MIRA with FASTA data
#  (corresponding qualities in .qual file loaded automatically)
#  -DIRECTORY:scf= optional, but SCFs can be found that way, so
#  that the automatic editor can work

# The sequences in the source were only masked by sequencing vector
#  but not subjected to quality clipping, MIRA can handle this
#  by doing a quality clip for itself

# we do NOT have an XML file with ancillary data, so also switch off
#  the option to load it

echo "Running mira"
mira -job=denovo,genome,normal,sanger -lowqualitydata --fasta=../data/fasta_set1/U13small_m.fasta -DIRECTORY:scf=../data/scf SANGER_SETTINGS -CLIPPING:qc=yes -LR:mxti=no | tee log_assembly.txt
echo 
echo "If all went well, the results are in 'mira_out*' files/directories."
echo "Consult log_assmbly.txt and/or stderr if errors occured."
echo
echo "(and please read the comments on top of the 'runme.sh' script)"

