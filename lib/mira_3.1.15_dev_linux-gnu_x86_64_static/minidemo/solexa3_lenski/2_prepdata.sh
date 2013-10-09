#!/bin/sh

######################################################################
#######
####### Prepare paired-end Solexa downloaded from NCBI
#######
######################################################################

# srrname:    is the SRR name as downloaded form NCBI SRA
# numreads:   maximum number of forward (and reverse) reads to take from
#              each file. Just to avoid bacterial projects with a coverage
#              of 200 or so.
# strainname: name of the strain which was re-sequenced

srrname="SRR030257"
numreads=5000000
strainname="REL8593A"

################################
# this is specific for this demo

mkdir data
cd data
ln -s ../origdata/NC_012967.gbk REL8593A-5000000_backbone_in.gbf

################################

numlines=$((4*${numreads}))

# put "/1" Solexa reads into file
echo "Copying ${numreads} reads from _1 (forward reads)"
zcat ../origdata/${srrname}_1.fastq.gz | head -${numlines} | sed -e 's/SRR[0-9.]*/&\/1/' >${strainname}-${numreads}_in.solexa.fastq

# put "/2" Solexa reads into file
echo "Copying ${numreads} reads from _2 (reverse reads)"
zcat ../origdata/${srrname}_2.fastq.gz | head -${numlines} | sed -e 's/SRR[0-9.]*/&\/2/' >>${strainname}-${numreads}_in.solexa.fastq

# make file with strainnames
echo "Creating file with strain names for copied reads (this may take a while)."
grep "@SRR" ${strainname}-${numreads}_in.solexa.fastq | cut -f 1 -d ' ' | sed -e 's/@//' -e "s/$/ ${strainname}/" >>${strainname}-${numreads}_straindata_in.txt

cd ..
