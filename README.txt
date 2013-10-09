###########################################################
# README: 
# NeatFreq version 1.0 - 10/1/2013
# J. McCorrison (ICS, J. Craig Venter Inst.)
# contact : jmccorri@jcvi.org
###########################################################

# 1.0 - INSTALL ##########################################

Deposit the entire contents of the NeatFreq distributable bundle in your desired install location. This bundle includes additional 1st and 3rd party tools available by open source license, expected in the NeatFreq install's /lib/ directory. 

(Optional) If you would like to remove the requirement to point to this area using the '-N' flag on every execution, update the global variable "$NEATFREQ_INSTALL" at the top of scripts NeatFreq_preprocess.pl, NeatFreq.pl, NeatFreq_postprocess_random.pl and NeatFreq_auto.pl.


# 2.0 - PREPROCESSING ####################################

NeatFreq assumes that only a single organism is present in the input dataset and the presence of a contaminant in any input data set may cause abnormal reduction of the target sample and maximal recruitment of the contaminant.

Targeted bin reduction with "-b align" requires preliminary read correction to remove as many false low frequency kmers as possible before NeatFreq preprocessing.  An example preprocessing pipeline is supplied in the assembly_automaton/ directory of the NeatFreq distributable code base.

High quantities of mate pairs may slow processing during targeted bin selection.  If this occurs, consider processing pairs as fragments, extracting as mate pairs following the process using NeatFreq_postprocess_random.pl.



# 3.1 - AUTOMATIC EXECUTION ##############################

Auto-run through each stage, automatically handling fasta/fastq file types and maintaining maximal paired end information.  To run a low memeory reduction of paired data without focused retention of 2-sided mates, enter all pairs as fragments.

EXAMPLE USAGE:
NeatFreq_auto.pl -N /install_loc/NeatFreq -z -q -k 19 -p PREFIX -m 1000 -x 10 -b align -f 33 -r fragments.fastq -rpairs pairs.fastq



# 3.2 - INDIVIDUAL PROCESSING EXECUTION ##################

##### 3.2.A - Standard Execution

EXAMPLE USAGE (preprocess):
NeatFreq_preprocess.pl -q -f 33 -m 19 -frag fragments.fastq -pair pairs.fastq

EXAMPLE USAGE (reduce):
NeatFreq.pl -p align_80 -b align -x 80 -r AUTO_FORMAT.frg.fasta -rpairs AUTO_FORMAT.prs.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m 10000 -v


##### 3.2.B - Using Pairs as Fragments

EXAMPLE USAGE (preprocess):
NeatFreq_preprocess.pl -q -f 33 -m 19 -pair pairs.fastq

EXAMPLE USAGE (reduce):
NeatFreq.pl -p align_as_fragments_80 -b align -x 80 -r AUTO_FORMAT.frg.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m 10000 -v

EXAMPLE USAGE (postprocess):
(Calculate number of fragments in AUTO_FORMAT.fasta for use in flag -numfrg):
NeatFreq_postprocess_random.pl -ids REDUCED_COVERAGE_ID_LIST.txt -numfrg 500 -frag AUTO_FORMAT.frg.fasta -pair AUTO_FORMAT.prs.fasta
