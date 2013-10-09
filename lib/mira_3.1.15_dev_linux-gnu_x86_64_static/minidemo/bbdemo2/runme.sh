#!/bin/sh

# this is a demo for a nice and simple genome MAPPING assembly with data which
# was downloaded from the NCBI trace archive, we have ancillary information in
# a XML file and TIGR made the reads, so we use that naming scheme

# using the "accurate" assembly qualifier turns on MIRA in high quality mode,
#  ready to tackle a fair number of problems one can encounter in genome
#  assembly.  
# It's not really needed for this small set, but neither does it harm.

# The data itself is pretty good, but has a few reads that show distinct signs
#  "extended too long", i.e., they contain really bad quality. This is why we
#  turn on -CL:bsqc with standard parameters (look at the parameters.par file)

# the difference to bbdemo1: we map reads from C.jejuni RM1221 against 
#  a backbone (the first 40kb of C.jejuni NCTC1168)

ln -f -s ../data/bbdataset1/cjejuni_demo* .

echo "Running mira"
mira -parameters=parameters.par | tee log_assembly.txt
echo "Done."
echo
echo -n "Using convert_project create files with info on SNPs betwee strains ... "
cd cjejuni_demo_assembly/cjejuni_demo_d_results
convert_project -f maf -t asnp cjejuni_demo_out.caf cjejuni_snpanalysis >/dev/null
cd ..
echo " done."
echo
echo
echo "Load the project into the GAP4 editor to have a look at the first 40kb of"
echo " Campylobacter jejuni RM1221 mapped against Campylobacter jejuni NCTC1168"
echo
echo "There are also showcase files to demonstrate how MIRA together with convert_project"
echo " helps you to find and analyse SNPs in prokaryotes."
echo "Normally, I would first load the project in gap4 (using caf2gap), perform a bit of"
echo " cleanup and manually check the data. Then only I would use gap2caf and"
echo " convert_project to extract and analyse SNP. But for showcasing, this will do."
echo
echo
echo "Want to know how many SNPs are between the strains?"
echo "  -> Look at the file cjejuni_snpanalysis_info_snplist.txt"
echo "Want to know the effect of SNPs on a gene/protein?"
echo "  -> Look at the file cjejuni_snpanalysis_info_featureanalysis.txt"
echo "Want to have a summary which genes/proteins are affected by SNPs?"
echo "  -> Look at the file cjejuni_snpanalysis_info_featuresummary.txt"
echo "Want to have the strain sequence of each gene/protein?" 
echo "  -> Look at the file cjejuni_snpanalysis_info_featuresequences.txt"
echo
echo "(and please read the comments on top of the 'runme.sh' script)"





#echo "HTML output was also switched on so that you can load 'cjejuni_demo_out.html'"
#echo " into a CSS compliant browser."
#echo ""
#echo "Look for SNPs between NCTC1168 (represented as 'read' named NC_002163) and"
#echo " reads from RM1221: search for the SROc tags which MIRA conveniently set for you."