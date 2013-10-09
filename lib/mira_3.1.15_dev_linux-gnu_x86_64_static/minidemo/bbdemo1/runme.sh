#!/bin/sh

# this is a demo for a nice and simple genome assembly with data which was
# downloaded from the NCBI trace archive, we have ancillary information in a
# XML file and TIGR made the reads, so we use that naming scheme

# using the "accurate" assembly qualifier turns on MIRA in high quality mode,
#  ready to tackle a fair number of problems one can encounter in genome
#  assembly.  
# It's not really needed for this small set, but neither does it harm.

# The data itself is pretty good, but has a few reads that show distinct signs
#  "extended too long", i.e., they contain really bad quality. This is why we
#  turn on -CL:bsqc with standard parameters



ln -f -s ../data/bbdataset1/cjejuni_demo* .

echo "Running mira"
mira -fasta -project=cjejuni_demo --job=denovo,genome,accurate,sanger SANGER_SETTINGS -LR:rns=TIGR -CL:bsqc=yes | tee log_assembly.txt
echo
echo
echo "This was a de-novo assembly of the first 40kb of Campylobacter jejuni RM1221."
echo "Load the project into the GAP4 editor to have a look at it."
echo
echo "(and please read the comments on top of the 'runme.sh' script)"




