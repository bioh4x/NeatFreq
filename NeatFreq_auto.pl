#!/usr/bin/perl -w
# J. Craig Venter Institute
# NeatFreq ver 1.0
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)

use strict;
use File::Basename; #required for fileparse()
use Getopt::Std; #better command line input method
use Getopt::Long; #better command line input method

my $NEATFREQ_INSTALL = "/home/jmccorri/scripts/NF_prod_testing";
my %Opts;

if (-e("./run.LOG.txt")){ system("rm ./run.LOG.txt"); }
system("touch ./run.LOG.txt");
open LOG, ">> ./run.LOG.txt";


################################################################################################################################################################################################################################################
# runsys : 
################################################################################################################################################################################################################################################
sub runsys {
	my $cmd = shift;
	print "\n\n++++ RUN CMD (NEATFREQ_AUTO) : \n++++ $cmd\n\n";
	print LOG "++++ RUN CMD (NEATFREQ_AUTO) : \n++++ $cmd\n\n";
	system("$cmd");
}

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


################################################################################################################################################################################################################################################
# outhelp : 
# print anything once to STDOUT and once to the LOG file -- pretty useful, eh?
################################################################################################################################################################################################################################################
sub outhelp {
	print "USAGE :\nNeatFreq_auto.pl (NeatFreq wrapper script)";
	printl "\n\t[ Count Mers -> Run Intelligent Selection -> Post-process to use maximal mate information ]\n\n";
	
	#SHARED
	print "REQUITED:\n";
	print "\t-p <prefix>\n";
	print "\t-r <reads.fasta>";
	print "\t-rpairs <interleaved.fasta/q> \n";
	print "\t\tFlags -r and -rpairs may be used alone for single library testing\n\t\tAll mates assumed to use inny orientation.\n";
	#PRE
	print "\t-b <bin selection method>\n";
	print "\t\t+ random = pull randomly from set (fast)\n\t\t+ align = run internal MSA within bins and include more low population sets (slow, -m flag required!)\n";
	print "\t\t-m <cd-hit-est memory cutoff in MB> (required for use with -b target)\n";
	print "\t-x <RMKFcutoff, set to approx. half of desired output coverage)\n";
	print "\t-q <toggle on fastq mode> (default = fasta)\n";
	print "\t-k <kmer size> (default = 19)\n";
	#NF
	print "\nOPTIONAL:\n";
	print "\t-v (silence print to screen, always verbose in logs)\n";
	print "\t-z = keep all output files (bin information), normally cleaned at the end of each run\n";
	print "\t-N = provide NeatFreq install location (if not updated in script)\n\n";
	print "\nExiting.";
	
	close LOG;
	exit;
}

################################################################################################################################################################################################################################################
sub bp_count {
	my $file = shift;
	my $inputcount = `/usr/local/packages/clc-ngs-cell/sequence_info $file | grep \'Total\' | awk \'{print \$2}\'`;
	chomp $inputcount;
	return $inputcount;
}

################################################################################################################################################################################################################################################
# MAIN : 
################################################################################################################################################################################################################################################
MAIN : {	
	my ($reads_file, $counts_file, $prefix, $cov_in, $bin_extract_type, $mem, $rpairs, $kmer_size, $offset, $NEW_NEATFREQ_INSTALL);
	my $pairstatus = 0;
	my $log = 1;
	my $keep = 0;
	my $fastq_bool = 0;
	my $status = GetOptions(\%Opts, "help!", "h!", "v!", "z!", "q!", 'k=s'=> \$kmer_size, 'N=s'=> \$NEW_NEATFREQ_INSTALL, 'f=s'=> \$offset, 'b=s'=> \$bin_extract_type, 'm=s'=> \$mem, 'r=s'=> \$reads_file, 'c=s'=> \$counts_file, 'p=s'=> \$prefix, 'x=s'=> \$cov_in, 'rpairs=s'=> \$rpairs);
	
	#PARSE INPUT------------------------------------------------------------------------------------------------------------------START
	print "\n++++ RUN CMD (NEATFREQ_AUTO) : NeatFreq Wrapper Started\n\n";
	
	if ( exists $Opts{help}){ printl("Help requested:\n"); outhelp(); }	
	if ( exists $Opts{h}){ printl("Help requested:\n"); outhelp(); }
	if ( exists $Opts{v}){ 
		printl("Print of logs to screen disabled.\n");
		$log = 0;
	}
	if ( exists $Opts{z} ){ 
		printl("Saving all output per user request (-z)\n");
		$keep = 1;
	}
	if (($NEW_NEATFREQ_INSTALL) && ($NEW_NEATFREQ_INSTALL ne $NEATFREQ_INSTALL)){
		printl("Forcing use of user-selected NeatFreq install directory : $NEW_NEATFREQ_INSTALL\n");
		if (!(-s("$NEW_NEATFREQ_INSTALL/NeatFreq.pl"))){
			printl("\n\nERROR! : No NeatFreq install found in suggested install directory.  \n\nExiting...\n");
			exit;
		}
		$NEATFREQ_INSTALL = $NEW_NEATFREQ_INSTALL;
	}
	if ( exists $Opts{q} ){ 
		printl("Fastq mode requested\n");
		$fastq_bool = 1;
	}
	if ( $kmer_size ){
		printl("Using User-Input Kmer Size : $kmer_size\n");
	}else{
		printl("Using Default Kmer Size : 19\n");
		$kmer_size = 19;
	}
	if ( $prefix ){
		printl("Using User-Input Prefix : $prefix\n");
	}else{
		printl("Using Default Prefix : OUT\n");
		$prefix = "OUT";
	}
	if ( $offset ){
		printl("Using User-Input Offset : $offset\n");
	}else{
		#no change
	}
	if (($cov_in) && ($cov_in > 0)){
		printl("Using input coverage : $cov_in\n");
	}else{
		printl("No input for -x <coverage> or input < 0.  EXITING!\n\n");
		outhelp();
	}
	if ($bin_extract_type){
		if (($bin_extract_type eq "random") || ($bin_extract_type eq "align")){
			printl("Using user-input bin extraction method = $bin_extract_type\n");
		}else{
			printl("ERROR:\nInput for flag -b does not match one of the allowed options (random or align).\nEXITING!\n\n");
			outhelp();
		}
	}else{
		printl("No input for -b <bin extract method>.\nUsing default value = random\n");
		$bin_extract_type = "random";
	}
	
	if ($mem){
		if ($bin_extract_type eq "align"){
			printl("Using memory cutoff for alignment within bins = $mem MB\n");
		}else{
			printl("WARNING: Memory flag will not be used because [-b align] was not used.\n");
		}
	}else{
		if ($bin_extract_type eq "align"){
			printl("ERROR: Memory flag not used with flag [-b align].\n");
			exit;
		}else{
			#ignore
		}
	}
	
	if ($reads_file && (-s("$reads_file")) ){
		printl("Using input reads file : $reads_file\n");
	}else{
		printl("WARNING:\tNo fragment-only read file input or reads files are 0 bytes in size. ( $reads_file )\nEXITING!\n");
		# outhelp();
	}
	
	if ( $rpairs && (-s("$rpairs"))){
		printl("Paired input reads input: $rpairs\n");	
		$pairstatus=1;
	}else{
		printl("WARNING:\tNo paired read file input or reads files are 0 bytes in size. ( $rpairs !\n");
		# outhelp();
	}
	
	my $input_status;
	# ver.P update : CHECK FOR READ STATUS
	if ( ( $rpairs && (-s("$rpairs"))) && ($reads_file && (-s("$reads_file"))) ){
		# both exist
		$input_status=2;
		printl("INPUT STATUS = 2, input sequence types found (pairs + fragments).\n");
	}elsif( (!(-s("$reads_file"))) ){
		# pairs only
		printl("INPUT STATUS = 1, Only PAIRED ENDs given as input.\n");
		$input_status=1;
	}elsif( (!(-s("$rpairs"))) ){
		# fragments only
		printl("INPUT STATUS = 0, Only FRAGMENTs given as input.\n");
		$input_status=0;
	}else{
		# neither - true fail case
		printl("No read file input or reads files are 0 bytes in size. ( $rpairs )\nEXITING!\n\n");
		# outhelp();
	}
	
		
	if ($NEW_NEATFREQ_INSTALL ne $NEATFREQ_INSTALL){
		#PREPROCESS
			if($input_status == 2){
				if ($fastq_bool == 1){ 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -q -f $offset -m $kmer_size -pair $rpairs -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -q -m $kmer_size -pair $rpairs -frag $reads_file"); }
				}
				else { 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -f $offset -m $kmer_size -pair $rpairs -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -m $kmer_size -pair $rpairs -frag $reads_file"); }
				}
			}elsif($input_status==1){
				if ($fastq_bool == 1){ 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -q -f $offset -m $kmer_size -pair $rpairs"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -q -m $kmer_size -pair $rpairs"); }
				}else{
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -f $offset -m $kmer_size -pair $rpairs"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -m $kmer_size -pair $rpairs"); }
				}
			}else{
				if ($fastq_bool == 1){ 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -q -f $offset -m $kmer_size -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -q -m $kmer_size -frag $reads_file"); }
				}
				else { 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -f $offset -m $kmer_size -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -m $kmer_size -frag $reads_file"); }
				}
			}
			
		# NEATFREQ EXECUTION
		
			# PREFIX
			my $NF_prefix = "$prefix" . "." . "$bin_extract_type" . "." . "$cov_in";
				
			if ((-s("AUTO_FORMAT.frg.fasta")) && (-s("AUTO_FORMAT.prs.fasta")) && (-s("allreads.mer_counts.minocc1.FIX.txt"))){
				runsys("$NEATFREQ_INSTALL/NeatFreq.pl -p $NF_prefix -b $bin_extract_type -x $cov_in -r AUTO_FORMAT.frg.fasta -rpairs AUTO_FORMAT.prs.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m $mem -v");
			}elsif(-s("AUTO_FORMAT.prs.fasta")) {
				runsys("$NEATFREQ_INSTALL/NeatFreq.pl -p $NF_prefix -b $bin_extract_type -x $cov_in -rpairs AUTO_FORMAT.prs.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m $mem -v");
			}elsif(-s("AUTO_FORMAT.frg.fasta")) {
				runsys("$NEATFREQ_INSTALL/NeatFreq.pl -p $NF_prefix -b $bin_extract_type -x $cov_in -r AUTO_FORMAT.frg.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m $mem -v");
			}else{
				printl("\n\nERROR! No output found from NeatfFreq_preprocess.pl - see NeatFreq_auto.log for details\n\n");
				exit;
			}
			
		#POSTPROCESS
			if ( $bin_extract_type eq "random" && $rpairs && (-s("$rpairs")) && $reads_file && (-s("$reads_file")) ){
					
				my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
				
				my $numfrg = `grep -c \'\>\' AUTO_FORMAT.frg.fasta`;
				if ( (-s("AUTO_FORMAT.frg.fasta")) && (!($numfrg > 0)) ){
					printl("\n\nERROR : Failure to parse fasta output of NeatFreq.pl.  See NeatFreq_auto.log for details.\n\n");
					exit;
				}else{
					$numfrg = 0;
				}
			
				runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -ids $ids -numfrg $numfrg -frag AUTO_FORMAT.frg.fasta -pair AUTO_FORMAT.prs.fasta");
			}
	}else{
		#PREPROCESS
			if($input_status == 2){
				if ($fastq_bool == 1){ 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -q -f $offset -m $kmer_size -pair $rpairs -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -q -m $kmer_size -pair $rpairs -frag $reads_file"); }
				}
				else { 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -f $offset -m $kmer_size -pair $rpairs -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -m $kmer_size -pair $rpairs -frag $reads_file"); }
				}
			}elsif($input_status==1){
				if ($fastq_bool == 1){ 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -q -f $offset -m $kmer_size -pair $rpairs"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -q -m $kmer_size -pair $rpairs"); }
				}else{
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -f $offset -m $kmer_size -pair $rpairs"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -m $kmer_size -pair $rpairs"); }
				}
			}else{
				if ($fastq_bool == 1){ 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -q -f $offset -m $kmer_size -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -q -m $kmer_size -frag $reads_file"); }
				}
				else { 
					if ($offset){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -f $offset -m $kmer_size -frag $reads_file"); }
					else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_preprocess.pl -N $NEW_NEATFREQ_INSTALL -m $kmer_size -frag $reads_file"); }
				}
			}
			
		# NEATFREQ EXECUTION
		
			# PREFIX
			my $NF_prefix = "$prefix" . "." . "$bin_extract_type" . "." . "$cov_in";
				
			if ((-s("AUTO_FORMAT.frg.fasta")) && (-s("AUTO_FORMAT.prs.fasta")) && (-s("allreads.mer_counts.minocc1.FIX.txt"))){
				runsys("$NEATFREQ_INSTALL/NeatFreq.pl -N $NEW_NEATFREQ_INSTALL -p $NF_prefix -b $bin_extract_type -x $cov_in -r AUTO_FORMAT.frg.fasta -rpairs AUTO_FORMAT.prs.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m $mem -v");
			}elsif(-s("AUTO_FORMAT.prs.fasta")) {
				runsys("$NEATFREQ_INSTALL/NeatFreq.pl -N $NEW_NEATFREQ_INSTALL -p $NF_prefix -b $bin_extract_type -x $cov_in -rpairs AUTO_FORMAT.prs.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m $mem -v");
			}elsif(-s("AUTO_FORMAT.frg.fasta")) {
				runsys("$NEATFREQ_INSTALL/NeatFreq.pl -N $NEW_NEATFREQ_INSTALL -p $NF_prefix -b $bin_extract_type -x $cov_in -r AUTO_FORMAT.frg.fasta -c allreads.mer_counts.minocc1.FIX.txt -z -m $mem -v");
			}else{
				printl("\n\nERROR! No output found from NeatfFreq_preprocess.pl - see NeatFreq_auto.log for details\n\n");
				exit;
			}
			
		#POSTPROCESS
			if ( $bin_extract_type eq "random" && $rpairs && (-s("$rpairs")) && $reads_file && (-s("$reads_file")) ){
					
				my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
				
				my $numfrg = `grep -c \'\>\' AUTO_FORMAT.frg.fasta`;
				if ( (-s("AUTO_FORMAT.frg.fasta")) && (!($numfrg > 0)) ){
					printl("\n\nERROR : Failure to parse fasta output of NeatFreq.pl.  See NeatFreq_auto.log for details.\n\n");
					exit;
				}else{
					$numfrg = 0;
				}
			
				runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -N $NEW_NEATFREQ_INSTALL -ids $ids -numfrg $numfrg -frag AUTO_FORMAT.frg.fasta -pair AUTO_FORMAT.prs.fasta");
			}
	}
		
	#CLEANUP
	printl("CLEANING UP!\n");
	#system("rm reads.des");
	#system("rm reads.sds");
	#system("rm reads.esq");
	#system("rm reads.ssp");
	#system("rm reads.ois");
	#system("rm reads.llv");
	#system("rm reads.suf");
	#system("rm reads.prj");
	#system("rm reads.lcp");
	#system("rm mock.fastq");
	#system("rm allreads.mer_counts.minocc1.txt");
	#system("rm tyr*");
	#system("rm extractFasta.log");
	#system("rm uniq_counts_ids.txt");
	#system("rm usable_reads_1.fasta");
	#system("rm AUTO_FORMAT*");
	system("rm -rf KMER_BINS_mp_frg KMER_BINS INITBINS_CDHITEST");
	#system("rm allreads.fastq");
	#system("rm allreads.fasta");
	#system("rm allreads.fasta.qual");
	#system("rm TESTreduced*");
	#system("mkdir run.logs");
	#system("mv *LOG* run.logs/");
	my $outprefix = "$prefix" . ".NeatFreq_auto.LOG.txt";
	system ("mv ./run.LOG.txt $outprefix");
	
	#printl("\nRUN COMPLETE!\nSee output file : $outfilename\n\nExiting.\n");
		
	close LOG;
}