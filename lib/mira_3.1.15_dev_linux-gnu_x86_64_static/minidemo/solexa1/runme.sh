
ln -s ../data/solexa_eco_art/ecoli_* .

mira --project=ecoli --job=mapping,genome,normal,solexa --fasta -AS:nop=1 -SB:bsn=ecoli_k12_mg1655 -SB:lsd=yes:bft=gbf SOLEXA_SETTINGS -AL:mrs=60:shme=5  | tee log_assembly.txt

echo
echo "Please read the README to know what is in the results."