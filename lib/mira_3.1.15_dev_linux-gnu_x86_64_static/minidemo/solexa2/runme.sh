
ln -s ../data/solexa_eco_art/ecoli_*fasta* .

mira --project=ecoli --job=denovo,genome,normal,solexa --fasta | tee log_assembly.txt

echo
echo "This was a de-novo assembly of Solexa reads."