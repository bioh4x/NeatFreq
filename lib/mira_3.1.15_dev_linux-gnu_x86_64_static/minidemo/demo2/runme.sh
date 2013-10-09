#!/bin/sh

# this is a quick demo on how to use MIRA from the command line
# 
# This is the same as 'demo1' (well run with all major default options), 
#  EXCEPT that our input files are not in this directory and we'll give
#  the output files/directories a special project name.
#
# We'll name the project "demo2" and tell the assembler where to find 
#  the file of filenames with the experiment files, the directory
#  with the experiment and the SCF files. Also, redirect the output
#  to run.log
#

echo "Running mira"
mira --project=demo2 -job=denovo,genome,normal,sanger -lowqualitydata -LR:eq=scf SANGER_SETTINGS -LR:ft=fofnexp:mxti=no  COMMON_SETTINGS -FILENAME:fofnexpin=../data/exp_set1/fofn -DIRECTORY:exp=../data/exp_set1:scf=../data/scf | tee log_assembly.txt
echo
echo "If all went well, the results are in 'demo2_results' directory."
echo " also have a look at files in 'demo2_info', e.g. at the contigstats file.."
echo "Consult log_assmbly.txt and/or stderr if errors occured."
echo
echo "(and please read the comments on top of the 'runme.sh' script)"

