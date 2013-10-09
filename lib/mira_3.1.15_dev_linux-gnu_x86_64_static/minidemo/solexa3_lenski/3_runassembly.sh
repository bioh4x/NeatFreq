#!/bin/bash

-mkdir assemblies
-mkdir assemblies/1sttest
cd assemblies/1sttest

lndir ../../data

echo "Starting mapping assembly, this can take a while (30-60 minutes)"
mira --project=REL8593A-5000000 --job=mapping,genome,accurate,solexa -SB:lsd=yes:bsn=ECO_B_REL606:bft=gbf >&log_assembly.txt

cd REL8593A-5000000_assembly/REL8593A-5000000_d_results

echo "Creating files (tables, HTML) with info on SNPs, please wait"
convert_project -f maf -t asnp -t hsnp REL8593A-5000000_out.maf REL8593A-5000000 >&/dev/null

cd ..
echo "All done, have a look at the directory: `pwd`"
