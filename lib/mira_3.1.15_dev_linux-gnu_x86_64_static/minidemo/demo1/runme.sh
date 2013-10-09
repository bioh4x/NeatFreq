#!/bin/sh

# this is a quick demo on how to use MIRA from the command line
# 
# easiest thing: all necessary files are in the directory where
#  we run mira and we use almost standard filenames and default options
#
# Note that this leads to a quite cluttered directory as is not
#  really recommended
#
# we'll make a denovo assembly, it's goinbg to be genomic data (no EST,
#  cDNA or transcripts), we want "normal" assembly and the sequences are
#  Sanger type
#
# The sequences are of low quality by todays standard, use -lowqualitydata
#
# just tell mira that 
#  - we load EXP files from a file of filenames and not the FASTA / XML 
#    combo that is used per default (-LR:sanft=fofnexp)
#  - switch off XML loading (-LR:mxti=no)
#  - load additional qualities from SCF as these EXP do not contain them
#    (-LR:eq=scf)
#
# The base qualities *are* important: leave out the eq=scf in the command line
#  and see how MIRA build then 2 instead of the correct 3 contigs

echo "Running mira"
mira -job=denovo,genome,normal,sanger -lowqualitydata -LR:eq=scf SANGER_SETTINGS -LR:ft=fofnexp:mxti=no | tee log_assembly.txt
echo
echo "If all went well, the results are in '*_results' directory."
echo " also have a look at files in '*_info', e.g. at the contigstats file.."
echo "Consult log_assmbly.txt and/or stderr if errors occured."
echo
echo "(and please read the comments on top of the 'runme.sh' script)"

