#!/bin/sh

# this is a quick demo on how to use MIRA with a parameter file

# This is the same as 'demo2' (well run with all major default options), 
#  EXCEPT that we use a parameter file to read the options of MIRA

echo "Running mira"
mira --params=demo3.par | tee log_assmbly.txt
echo
echo "If all went well, the results are in 'demo3_out*' files/directories."
echo "Consult log_assmbly.txt and/or stderr if errors occured."
echo
echo "(and please read the comments on top of the 'runme.sh' script)"


