###########################################################
# README: 
# Assembly Automaton version 1.0 - 10/1/2013
# J. McCorrison (ICS, J. Craig Venter Inst.)
# contact : jmccorri@jcvi.org
###########################################################

# 1.O - Install ###########################################

This software is dependent on a number of 3rd party utilities and requires updated to the hardcoded install locations at the top of the script (assembly_automaton.pl).

3rd party tools requiring their own install and configuration include:
* MIRA
* CLC NGS Cell
* Cutadapt
* Quake
* Python
* Perl
* Fastq Utilities
* Dust
* Rnnotator
* Celera WGS Assembler
* Velvet Assembler
* Newbler Assembler



# 2.O - How It Works ######################################

The "Assembly Automaton" was originally built for high throughput pre-processing of data containing extremely deep and/or highly variable coverage.  Users may pick and choose the processing methods they require (always run in the order of the pipeline depicted below).  The tool will track all sequencing data as it passes through each processing step and log information on the quantity and paired end percentages of reads after each step.  A brief CLC assembly step is used following every process, allowing the user to evaluate the effect each process had on the resulting consensus calls.  Following all pre-processing, multiple assemblies are run and statistics commands are delivered for a user to run for evaluation of the assembled result.

Serial checks are as follows, using the tools described below each processing stage:

    Contaminant Checks
        clc_ref_assemble_long to contaminant references using 40% query sequence length and 95% identity minimum match requirements
    Read Correction (illumina-only)
        quake.py – very long run times on high coverage samples (or)
        allpaths read correction module – corrects fragments as fake mates
    Unique Check
        fastq::dedupe - not compatible with allpaths!
    QV Trim, Rapid Assemblies, Read Alignment, Sequence Info :
        CLC NGS cell
    Assembly :
        clc_novo_assemble, Velvet, CA, Newbler
    Low Complexity Check :
        DUST
    Adapter Removal :
        Cutadapt – run in serial on input adapters
            Forward and reverse orientation with random hexamers added
            Front and Back Match (run in serial) : Trim at as few as 6 bp, Use output > 50 bp
            Middle Match (run last) : throw out matches (likely chimeric)
    Data conversions :
        gatekeeper, sffToCA, fastqToCA
		

Notes and Limitations of version 1.0:

    Illumina only processes : read correction, quality trim
    Run regardless of flag inclusion : quality trim (illumina-only)
    Not yet included : SPAdes, IBDA-UD assembly
    Runs in serial on 1 host (initial build is for small bacterial genomes)
    To use pre-processing output without assembly, you must kill the automaton job once velvet assembly has started.
    This script is still in development - please contact support with any questions, errors or requests!

Supported Inputs :

    SFF fragments
    SFF pairs
    Illumina fragments
    Illumina paired ends but NOT mate pairs (orientation currently hardcoded)



# 3.O - Example Usage ###################################

## Automaton #######

$NEATFREQ_INSTALL/assembly_automaton/lib/shuffleSequences_fastq.pl seq.1.txt seq.2.txt interleaved_mates.fastq

$NEATFREQ_INSTALL/assembly_automaton/assembly_automaton.pl -p PREFIX -g auto -d -r allpaths -f 33 -o /output/directory/ -precontam contaminant1.fasta -adapter GCCGGAGCTCTGCAGATATC -adapter CGAGAGATACTGTACTAGAGCG -contam contaminant2.fasta -SFFfragment fragments.sff -ILLUMfragment fragments.fastq -ILLUMpair interleaved_mates.fastq -ILLUMinsert 260


## Output Files ####

Output files from processing will be supplied in:

    FASTQ = file locations listed in $output_directory/$PREFIX.post_automaton_fastq_locs.txt
    FASTA (Newbler formatted) - $output_directory/NEWBLER_ASSEMBLIES/*fna
    FRG (Celera formatted) - /usr/local/depot/projects/HMP1/HMPMDA_BULK.tmpdir/HMPMDA0240/AUTOMATON/CELERA_ASSEMBLIES/*frg


## Read Tracking ####

Read tracking at each stage can be found in the output logs within $cwd/assembly_automaton.LOG.txt

Quickly extract sequence metadata following each processing stage using grep : ( grep -A 5 'Perc' $cwd/assembly_automaton.LOG.txt )

