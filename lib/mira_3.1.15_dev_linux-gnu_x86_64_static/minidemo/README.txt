IMPORTANT NOTICE 1: if you're looking for 454 demos, please follow the
walkthrough guides in the documentation on 454 assembly



IMPORTANT NOTICE 2: (only applicable to source code distributions). Please do
NOT run mira in the original 'minidemo' directory if you plan to make a binary
("make distrib") or source distribution ("make dist") of MIRA! Doing so would
mean that result and temporary files would get included in the resulting
package.



This directory contains a few demos on how to use the MIRA assembler. The
different demos are all in sub-directories called 'demo1', 'demo2',
'estdemo1', 'bbdemo1' etc. and demonstrates basic usage.

Make sure that you have the binaries of the MIRA distribution in your
path. Then simply change to any demo directory and start the
'runme.sh' scripts present there. Have a look at them to see which
parameters they use, have also a look at the parameter files (if
present) in the directories.

Once an assembly is finished, you will find three directories in the demo
directory where you started (where "..." stands for a variable prefix):

1) ..._results: this directory contains all the output files of the assembly
   in different formats. there are many different formats that can be chosen
   from the command line: fasta, caf, ace, gap4 directed assembly, simple text
   etc.
   The following files in ..._results contain results of the assembly in
   different formats:

   a)..._out.txt: this file contains in a human readable format the aligned
    assembly results, where all input sequences are shown in the context of the
    contig they were assembled into. This file is just meant as a quick way for
    people to have a look at their assembly without specialised alignment
    finishing tools.
   b)..._out.padded.fasta: this file contains as FASTA sequence the consensus
    of the contigs that were assembled in the process. Positions in the
    consensus containing gaps (also called 'pads', denoted by an asterisk) are
    still present. The computed consensus qualities are in the corresponding
    ..._out.padded.fasta.qual file.
   c)..._out.unpadded.fasta: as above, this file contains as FASTA sequence
    the consensus of the contigs that were assembled in the process, put
    positions in the consensus containing gaps were removed. The computed
    consensus qualities are in the corresponding ..._out.unpadded.fasta.qual
    file.
   d) ..._out.caf: this is the result of the assembly in CAF format, which can
    be further worked on with, e.g., tools from the caftools package from the
    Sanger Centre and later on be imported into, e.g., the Staden gap4 assembly
    and finishing tool.
   e) ..._out.ace: this is the result of the assembly in ACE format. This
    format can be read by viewers like the TIGR clview or by consed from the
    phred/phrap/consed package.
   f)..._out.gap4da: this directory contains the result of the assembly suited
    for the direct assembly import of the Staden gap4 assembly viewer and
    finishing tool.


2) ..._info: this directory contains information files of the final
   assembly. They provide statistics as well as, e.g., information (easily
   parseable by scripts) on which read is found in which contig etc. 

   The following files in bchoc_info contain statistics and other results of
   the assembly:

   a) ..._info_callparameters.txt: This file contains the parameters as given
    on the mira command line when the assembly was started.
   b) ..._info_contigstats.txt: This file contains in tabular format
    statistics about the contigs themselves, their length, average consensus
    quality, number of reads, maximum and average coverage, average read
    length, number of A, C, G, T, N, X and gaps in consensus.
   c) ..._info_contigreadlist.txt: This file contains information which reads
    have been assembled into which contigs (or singlets).
   d) ..._info_consensustaglist.txt: This file contains information about the
    tags (and their position) that are present in the consensus of a contig.
   e) ..._info_readstooshort: A list containing the names of those reads that
    have been sorted out of the assembly only due to the fact that they were
    too short, before any processing started.
   f) ..._info_readtaglist.txt: This file contains information about the tags
    and their position that are present in each read. The read positions are
    given relative to the forward direction of the sequence (i.e. as it was
    entered into the the assembly).
   g) ..._error_reads_invalid: A list of sequences that have been found to be
    invalid due to various reasons (given in the output of the assembler).

3) ..._log: this directory contains log files and temporary assembly files. It
   can be safely removed after an assembly as there may be easily a few GB of
   data in there that are not normally not needed anymore.


Example of repeat resolving and different ways to call MIRA
-----------------------------------------------------------

It's a toy project, just to show a few things

Located in the directories demo1 to demo5. Change to any demo directory and
start the 'runme.sh' scripts present there.

The dataset used for these demos is contained in the directories
'data/exp_set1' (for a Staden experiment file set), 'data/scf' (for
all SCF trace files), and 'data/fasta_set1' (for the same dataset
converted to FASTA format).

The dataset itself is extremely small (26 sequences only) but contains some
tricky problems, like inverse repeats and SNPs due to samples sequenced from
different organisms (different humans in this case). The data was generated
some years ago on ABI373 sequencer machines and is taken from a larger
project. Only the regions with the repeats are included, the sequences do NOT
form a single contig. Depending on the dataset and how it is treated (quality
and/or SCFs available, clipping etc.) and the assembler (phrap, cap2/3/4, pga,
mira), one gets 2, 3 or 4 contigs. 2 contigs is plain wrong, different repeats
are then assembled together (most assemblers will give you that). 3 contigs is
correct (MIRA 2.8.x did that). Unfortunately, MIRA 2.9.x and later make 3
contigs and one singlet out of the data set. This is due to more extended
repeat recognition routines: a base position that is crucial for
disambiguating repeats is also containing sequencing errors in two of the
reads ... which then also get, due to sequencing errors, wrongly used for
contig disambiguation.

Some MIRA parameter settings in the demos produce intentionally wrong
assemblies (e.g. when repeat marking is switched off) to give an idea
on what some options actually do.



Larger genome assembly and mapping assembly:
--------------------------------------------

Two examples for genome assembly. One for genome assembly from scratch
(bbdemo1) and one with genome assembly against a backbone - a MAPPING assembly
- in bbdemo2. The later is especially interesting as it shows how to perform
mapping assemblies and how MIRA recognises SNPs between strains while
differentiating them from SNPs within a strain, but in different repeats.

The dataset used for these demos is contained in the directories
'data/bbdataset1'. These sequences and ancillary information for Campylobacter
jejuni RM1221 contained in these files were downloaded from the NCBI trace
archive and put together into one FASTA (and .qual) file.

For the backbone in 'bbdemo2', the first 40 kb of the genome of Campylobacter
jejuni NCTC 11168 were downloaded from the NCBI as GenBank file (NC_002163.1).



De-novo assembly and mapping of Solexa data:
--------------------------------------------

Two toy examples with an artificial data set, just as demo. solexa1 contains a
mapping assembly, solexa2 a de-novo. Especially the first should be
interesting as it shows how to perform mapping assemblies and how MIRA
recognises SNPs between strains.

One additional example: solexa3_lenski

That demo allows you to perform a mapping assembly with data from Richard
Lenski deposited at the NCBI. For more info, consult the walkthrough in the
MIRA help manual on Solexa data.



EST assembly and SNP analysis:
------------------------------

The data is public, taken from NCBI (see README in data directory) and
demonstrates how miraEST assembles ESTs to transcripts and how
different SNP types are set.

You might want to start with estdemo2 as this is a really tiny
showcase with 6 sequences.

Then proceed to the first demo and assembly a small real life
project. Note that the preprocessing cleaning was not done optimally,
I do not have the real vector sequences that were used in that project.

If you have some bandwidth at your disposal, you also could download
the trace data from NCBI ... which I did not include for obvious
reasons as this package is already large enough without those, thank
you.

Other interesting things to try: looking at the differences of the
assemblies when traces or even quality files are left out.
