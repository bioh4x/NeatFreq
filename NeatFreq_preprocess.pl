#!/usr/bin/perl -w
#INITIAL SETUP
# J. Craig Venter Institute
# NeatFreq ver 1.0.1
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)

use strict;
use File::Basename; #required for fileparse()
use Getopt::Std; #better command line input method
use Getopt::Long; #better command line input method

#OPEN LOG
my $NEATFREQ_INSTALL = `pwd`;
chomp $NEATFREQ_INSTALL;
open LOG, "> ./NeatFreq_preprocess.LOG.txt";

########################################################################################################################
# printl : 
# print anything once to STDOUT and once to the LOG file -- pretty useful, eh?
########################################################################################################################
sub printl {
	my $pr = shift;
	#run STDOUT print
	print $pr;
	#run LOG print (when appropriate)
	print LOG $pr;
}


########################################################################################################################
# outhelp: 
# print help text
########################################################################################################################
sub outhelp {
	printl("NeatFreq_preprocess.pl\n");
	printl("\tRequired Flags:\n");
	printl("\t-frag fragments.fast[a/q]\n\tAND/OR\n\t-pair pairs.interleaved.fast[a/q]\n\n");
	printl("\tOptional Flags:\n");
	printl("\t-q (toggles fastq mode, must be used with fastq files!)\n");
	printl("\t-f <qv offset>\n");
	printl("\t-N = provide NeatFreq install location (if not updated in script)\n\n");
	
	exit;
}

################################################################################################################################################################################################################################################
# runsys : 
################################################################################################################################################################################################################################################
sub runsys {
	my $cmd = shift;
	printl("++ RUN CMD : $cmd\n");
	system("$cmd");
}


########################################################################################################################
# fixfasta: 
########################################################################################################################
sub fixfasta {
	my $fragfile = shift;
	my $pairfile = shift;
	my $fastq_status = shift;
	my $mersize = shift;
	
	printl("\n= 01 FIX IDS =\n");
	if( (-s($fragfile)) && (-s($pairfile)) ){
		runsys("$NEATFREQ_INSTALL/configure_automaton_for_neatfreq.pl -N $NEATFREQ_INSTALL -frg $fragfile -prs $pairfile -out AUTO_FORMAT");
		runsys("cat AUTO_FORMAT.frg.fasta > allreads.fasta");
		runsys("cat AUTO_FORMAT.prs.fasta >> allreads.fasta");
	}
	elsif(-s($fragfile)){
		runsys("$NEATFREQ_INSTALL/lib/fastx_renamer -n COUNT -i $fragfile -o allreads.fasta");
	}else{
		runsys("$NEATFREQ_INSTALL/lib/fastx_renamer -n COUNT -i $pairfile -o allreads.fasta");
	}
	
	printl("\n\n= 02 TALLYMER =\n");
	runsys("$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt suffixerator -dna -pl -tis -suf -lcp -lossless -v -parts 4 -db allreads.fasta -indexname reads");
	runsys("$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt tallymer mkindex -mersize $mersize -minocc 1 -indexname tyr-reads-minocc1 -counts -pl -esa reads");
	runsys("$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt tallymer search -strand fp -output qseqnum qpos counts sequence -tyr tyr-reads-minocc1 -q allreads.fasta > allreads.mer_counts.minocc1.txt");
	printl("\n\n= 03 CLEAN NONEXISTENT MERS (due to non-basecalls) =\n");
	runsys("$NEATFREQ_INSTALL/clean_up_zero_mers_4_mernalysis.pl allreads.mer_counts.minocc1.txt allreads.mer_counts.minocc1.FIX.txt");
}



########################################################################################################################
# fixfastq: 
########################################################################################################################
sub fixfastq {
	my $fragfile = shift;
	my $pairfile = shift;
	my $fastq_status = shift;
	my $mersize = shift;
	my $offset = shift;
	
	printl("\n= 01 FIX IDS AND CREATE FASTA EQUIVALENTS =\n");
	if( (-s($fragfile)) && (-s($pairfile)) ){
		if ($offset != 0){
			runsys("$NEATFREQ_INSTALL/configure_automaton_fastq_for_neatfreq.pl -N $NEATFREQ_INSTALL -f $offset -frg $fragfile -prs $pairfile -out AUTO_FORMAT");
		}else{
			runsys("$NEATFREQ_INSTALL/configure_automaton_fastq_for_neatfreq.pl -N $NEATFREQ_INSTALL -frg $fragfile -prs $pairfile -out AUTO_FORMAT");
		}
		runsys("cat AUTO_FORMAT.frg.fastq > allreads.fastq");
		runsys("cat AUTO_FORMAT.prs.fastq >> allreads.fastq");
		
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual allreads.fastq allreads.fasta allreads.fasta.qual");
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual AUTO_FORMAT.frg.fastq AUTO_FORMAT.frg.fasta AUTO_FORMAT.frg.fasta.qual");
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual AUTO_FORMAT.prs.fastq AUTO_FORMAT.prs.fasta AUTO_FORMAT.prs.fasta.qual");
		
	}elsif(-s($fragfile)){
		if ($offset != 0){
			runsys("$NEATFREQ_INSTALL/lib/fastx_renamer -Q $offset -n COUNT -i $fragfile -o allreads.fastq");
			if (!(-s("./allreads.fastq"))){ printl("ERROR : Offset incorrect (see error above, fix with -f flag).\n\nExiting..."); exit; }
		}else{
			runsys("$NEATFREQ_INSTALL/lib/fastx_renamer -n COUNT -i $fragfile -o allreads.fastq");
		}
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual allreads.fastq AUTO_FORMAT.frg.fasta AUTO_FORMAT.frg.fasta.qual");
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual allreads.fastq allreads.fasta allreads.fasta.qual");
	}else{
		if ($offset != 0){
			runsys("$NEATFREQ_INSTALL/lib/fastx_renamer -Q $offset -n COUNT -i $pairfile -o allreads.fastq");
			if (!(-s("./allreads.fastq"))){ printl("ERROR : Offset incorrect (see error above, fix with -f flag).\n\nExiting..."); exit; }
		}else{
			runsys("$NEATFREQ_INSTALL/lib/fastx_renamer -n COUNT -i $pairfile -o allreads.fastq");
		}
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual allreads.fastq AUTO_FORMAT.prs.fasta AUTO_FORMAT.prs.fasta.qual");
		runsys("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual allreads.fastq allreads.fasta allreads.fasta.qual");
	}
	printl("\n\n= 02 TALLYMER =\n");
	runsys("$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt suffixerator -dna -pl -tis -suf -lcp -lossless -v -parts 4 -db allreads.fasta -indexname reads");
	runsys("$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt tallymer mkindex -mersize $mersize -minocc 1 -indexname tyr-reads-minocc1 -counts -pl -esa reads");
	runsys("$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt tallymer search -strand fp -output qseqnum qpos counts sequence -tyr tyr-reads-minocc1 -q allreads.fasta > allreads.mer_counts.minocc1.txt");
	printl("\n\n= 03 CLEAN NONEXISTENT MERS (due to non-basecalls) =\n");
	runsys("$NEATFREQ_INSTALL/clean_up_zero_mers_4_mernalysis.pl allreads.mer_counts.minocc1.txt allreads.mer_counts.minocc1.FIX.txt");
}



########################################################################################################################
#       MAIN                                                                                                           #
########################################################################################################################
MAIN : {
#command line input vars
	my %Opts;
	my ($fragfile, $pairfile, $fastq_status, $mersize, $offset, $cleanup);
	my $NEW_NEATFREQ_INSTALL = `pwd`;
#pull in command line input
	my $status = GetOptions(\%Opts, "help", "h", "q", "z", 'm=s'=> \$mersize, 'f=s'=> \$offset, 'N=s'=> \$NEW_NEATFREQ_INSTALL, 'frag=s'=> \$fragfile, 'pair=s'=> \$pairfile);
	if ( exists $Opts{help} || exists $Opts{h}){ outhelp(); }
	if ( exists $Opts{q} ){ $fastq_status = 1; }
	else{ $fastq_status = 0; }
	if ($mersize){
		printl("\nUSING INPUT MER SIZE = $mersize\n");
	}else{
		printl("\nUSING DEFAULT MER SIZE = 19\n");
		$mersize=19;
	}
	if ($offset){
		printl("\nUSING INPUT QV OFFSET = $offset\n");
	}else{
		$offset = 0;
	}
	if (($NEW_NEATFREQ_INSTALL) && ($NEW_NEATFREQ_INSTALL ne $NEATFREQ_INSTALL)){
		printl("Forcing use of user-selected NeatFreq install directory : $NEW_NEATFREQ_INSTALL\n");
		if (!(-s("$NEW_NEATFREQ_INSTALL/NeatFreq.pl"))){
			printl("\n\nERROR! : No NeatFreq install found in suggested install directory.  \n\nExiting...\n");
			exit;
		}
		$NEATFREQ_INSTALL = $NEW_NEATFREQ_INSTALL;
	}
		if (!(-s("$NEW_NEATFREQ_INSTALL/NeatFreq.pl"))){
			printl("\n\nERROR! : No NeatFreq install found in suggested install directory. ( $NEW_NEATFREQ_INSTALL/NeatFreq.pl )\n\nExiting...\n");
			exit;
		}
	
	if ( ($fragfile) && ($pairfile) && (-s("$fragfile")) && (-s("$pairfile")) ){
		printl("\nUSING INPUT FRAGMENT-ONLY SEQUENCE FILE : $fragfile\n"); 
		printl("\nUSING INPUT PAIRED END SEQUENCE FILE : $pairfile\n"); 
		
		if ($fastq_status == 0){
			fixfasta($fragfile, $pairfile, $fastq_status, $mersize);
		}else{
			fixfastq($fragfile, $pairfile, $fastq_status, $mersize, $offset);
		}
		
	}
	elsif ( $fragfile ){ 
		if (-s ("$fragfile")) { 
			printl("\nUSING INPUT FRAGMENT-ONLY SEQUENCE FILE : $fragfile\n"); 
		}
		else{ 
			printl("\nPROBLEM FOUND : File $fragfile does not exist or has 0 byte size.\nExiting...");
			#exit;
		}  
		
		
		if ($fastq_status == 0){
			fixfasta($fragfile, "NONE", $fastq_status, $mersize);
		}else{
			fixfastq($fragfile, "NONE", $fastq_status, $mersize, $offset);
		}
		
	}
	elsif ( $pairfile ){ 
		if (-s ("$pairfile")) { 
			printl("\nUSING INPUT PAIRED END SEQUENCE FILE : $pairfile\n"); 
		}
		else{ 
			printl("\nPROBLEM FOUND : File $pairfile does not exist or has 0 byte size.\nExiting...");
			#exit;
		}  
		
		if ($fastq_status == 0){
			fixfasta("NONE", $pairfile, $fastq_status, $mersize);
		}else{
			fixfastq("NONE", $pairfile, $fastq_status, $mersize, $offset);
		}
	}else{
		printl("\nERROR : No input sequences found.  See usage:\n\n");
		outhelp();
	}
	
# CLEAN UP

	runsys("rm reads.des  reads.esq  reads.lcp  reads.llv  reads.ois  reads.prj  reads.sds  reads.ssp  reads.suf");
	runsys("rm tyr-reads-minocc1.mbd  tyr-reads-minocc1.mct  tyr-reads-minocc1.mer");
	runsys("rm allreads.mer_counts.minocc1.txt");
	
#Fix 'em
	
	if (-s("allreads.mer_counts.minocc1.FIX.txt")){
		if ($fragfile && $pairfile){
			printl("\n\nSuccess! See NeatFreq input...\n");
			printl("\t-r AUTO_FORMAT.frg.fasta\n");
			printl("\t-rpairs AUTO_FORMAT.prs.fasta\n");
		}
		elsif ($fragfile){ 
			runsys("mv allreads.fasta AUTO_FORMAT.frg.fasta");
			runsys("mv allreads.fastq AUTO_FORMAT.frg.fastq");
			printl("\n\nSuccess! See NeatFreq input...\n");
			printl("\t-r AUTO_FORMAT.frg.fasta\n"); 
		}
		elsif ($pairfile ){ 
			runsys("mv allreads.fasta AUTO_FORMAT.prs.fasta");
			runsys("mv allreads.fastq AUTO_FORMAT.prs.fastq");
			printl("\n\nSuccess! See NeatFreq input...\n");
			printl("\t-rpairs AUTO_FORMAT.prs.fasta\n"); 
		}
		printl("\t-c allreads.mer_counts.minocc1.FIX.txt\n\n");
		if ($fastq_status == 1){
			#printl("\nTo extract quality scores following NeatFreq, use commands:\n");
			#if ($fragfile){ printl("\t/usr/local/devel/BCIS/assembly/tools/extractFasta -i AUTO_FORMAT.frg.fastq -idlist DDD.REDUCED_COVERAGE_ID_LIST.txt -o NeatFreq_reduced_fragments.fastq\n"); }
			#if ($pairfile ){ printl("\t/usr/local/devel/BCIS/assembly/tools/extractFasta -i AUTO_FORMAT.prs.fastq -idlist DDD.REDUCED_COVERAGE_ID_LIST.txt -o NeatFreq_reduced_pairs.fastq\n"); }
		}
	}
	

	
	
	exit;
}
