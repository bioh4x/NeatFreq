#!/bin/sh

# this is essentially the same as demo4 EXCEPT that 
#  it works with completely unmasked sequences (i.e. no sequencing vector 
#  was masked, the complete reads are present) and demonstrates the
#  effectiveness of a number of clipping options 
#
# Though please note that this is not the recommended working mode of
#  any assembler (including MIRA), one should _really_ tag/mask/whatever
#  sequencing vector with specialised programs (e.g. pregap4, cross_match,
#  scylla and others).

echo "Running mira on non vector masked data in FASTA standard mode,"
echo " i.e. with all sorts of security clippings. The quality clips"
echo " must be defined extra. Please wait ..."
mira  --project=nmpvc -job=denovo,genome,normal,sanger -lowqualitydata --fasta=../data/fasta_set1/U13small_nm.fasta -DIRECTORY:scf=../data/scf SANGER_SETTINGS -CLIPPING:qc=yes -LR:mxti=no > run_nmpvc.log
echo 
echo "If all went well, the results are in 'nmpvc_results' directory."
echo "Consult run_nmpvc.log and/or stderr if errors occured."

echo 
echo
echo "Running mira on non vector masked data but only with quality,"
echo " clippings. Please wait ..."
mira --project=nmqconly --job=denovo,genome,normal,sanger -lowqualitydata --fasta=../data/fasta_set1/U13small_nm.fasta -DIRECTORY:scf=../data/scf SANGER_SETTINGS -CLIPPING:pvlc=no:qc=yes:emlc=no -LR:mxti=no > run_nmqconly.log
echo 
echo "If all went well, the results are in 'nmqconly_results' directory."
echo "Consult run_nmqconly.log and/or stderr if errors occured."

echo 
echo
echo "Running mira on non vector masked data but without any clipping"
echo " clippings. Please wait ..."
mira --project=nmnocl --job=denovo,genome,normal,sanger -lowqualitydata --noclipping --fasta=../data/fasta_set1/U13small_nm.fasta -DIRECTORY:scf=../data/scf -LR:mxti=no > run_nmnocl.log
echo 
echo "If all went well, the results are in 'nmnocl_results' directory."
echo "Consult run_nmnocl.log and/or stderr if errors occured."
echo
echo "(and please read the comments on top of the 'runme.sh' script)"
