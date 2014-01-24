#!/usr/bin/perl -w
# J. Craig Venter Institute
# NeatFreq ver 1.0
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)

use strict;
use File::Basename; #required for fileparse()
use Getopt::Std; #better command line input method
use Getopt::Long; #better command line input method

my $NEATFREQ_INSTALL = "/usr/local/devel/BCIS/assembly/tools/NeatFreq";
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
	print "\nUSAGE :\nNeatFreq_auto.pl (NeatFreq wrapper script)";
	print "\n\t[ Count Mers -> Run Intelligent Selection -> Post-process to use maximal mate information ]\n\n";
	
	#SHARED
	print "REQUIRED:\n";
	print "\t-p <prefix>\n";
	print "\t-r <reads.fasta>";
	print "\t-rpairs <interleaved.fasta/q> \n";
	print "\t\tFlags -r and -rpairs may be used alone for single library testing\n\t\tAll mates assumed to use inny orientation.\n";
	
	#PRE
	print "\t-b <bin selection method>\n";
	print "\t\t+ random = pull randomly from set (fast)\n\t\t+ align = run internal MSA within bins and include more low population sets (slow, -m flag required!)\n";
	print "\t\t-m <cd-hit-est memory cutoff in MB> (required for use with -b target)\n";
	print "\t-x <RMKFcutoff, set to approx. half of desired output coverage)\n";
	print "\t-g <run analysis as fragment-only (low memory) but extract mates from reads input as pairs>\n";
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
	my $run_as_frag_bool = 0;
	my $status = GetOptions(\%Opts, "help!", "h!", "v!", "z!", "q!", 'g!', 'k=s'=> \$kmer_size, 'N=s'=> \$NEW_NEATFREQ_INSTALL, 'f=s'=> \$offset, 'b=s'=> \$bin_extract_type, 'm=s'=> \$mem, 'r=s'=> \$reads_file, 'c=s'=> \$counts_file, 'p=s'=> \$prefix, 'x=s'=> \$cov_in, 'rpairs=s'=> \$rpairs);
	
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
	if ( exists $Opts{g} ){ 
		printl("Run all analysis as fragments (low memory) but recruit mates from input pairs.\n");
		$run_as_frag_bool = 1;
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
		printl("WARNING:\tNo fragment-only read file input or reads files are 0 bytes in size. ( $reads_file )!\n");
		# outhelp();
	}
	
	if ( $rpairs && (-s("$rpairs"))){
		printl("Paired input reads input: $rpairs\n");	
		$pairstatus=1;
	}else{
		printl("WARNING:\tNo paired read file input or reads files are 0 bytes in size. ( $rpairs !\n");
		# outhelp();
	}
	
	
	# CHECK FOR RUN AS FRAGMENTS
	my $orig_fragments;
	my $toggle = 0;
	if ($run_as_frag_bool == 1) {
		
		if ($fastq_bool == 1){ 
			if ($reads_file) {	
				system("$NEATFREQ_INSTALL/lib/fastq_to_fasta_qual $reads_file tmp.ALL_FRAGS.fasta tmp.ALL_FRAGS.qual"); 
				$orig_fragments = "tmp.ALL_FRAGS.fasta";
				if (!(-s("tmp.ALL_FRAGS.fasta"))){ printl("++EXITING. CONVERSION FAILED.\n"); exit; }
			}
		}else{
			if ($reads_file) {	$orig_fragments = $reads_file; }
		}
		
		if (( $rpairs && (-s("$rpairs"))) && ($reads_file && (-s("$reads_file"))) ){
			printl("+ WARNING : Running with fragments and pairs as fragment-only.  Merging sequences...\n");
			my $newrpairs = "./ALL_INPUT.";
			if ($fastq_bool == 1){$newrpairs = "$newrpairs" . "fastq";}
			else{$newrpairs = "$newrpairs" . "fasta";}
			runsys("cat $reads_file > $newrpairs");
			runsys("cat $rpairs >> $newrpairs");
			$toggle = 1;
			$reads_file = $newrpairs;
			undef($rpairs);
		}elsif ( $reads_file && (-s("$reads_file")) ){
			#no change required
		}elsif ( $rpairs && (-s("$rpairs")) ){
			$reads_file = $rpairs;
			undef($rpairs);
			$toggle = 1;
		}else{
			printl("No input sequences.\nExiting...\n");
			exit;
		}
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
				printl("\n\nERROR! No output found from NeatFreq_preprocess.pl - see NeatFreq_auto.log for details\n\n");
				exit;
			}
			
		#POSTPROCESS
			if ( $run_as_frag_bool == 0 && ($bin_extract_type eq "random") && $rpairs && (-s("$rpairs")) && $reads_file && (-s("$reads_file")) ){
					
				my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
				
				my $numfrg = `grep -c \'\>\' AUTO_FORMAT.frg.fasta`;
				chomp $numfrg;
				if ( (-s("AUTO_FORMAT.frg.fasta")) && (!($numfrg > 0)) ){
					printl("\n\nERROR : Failure to parse fasta output of NeatFreq.pl.  See NeatFreq_auto.log for details.\n\n");
					exit;
				}else{
					$numfrg = 0;
				}
				
				if ($fastq_bool == 0){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -ids $ids -numfrg $numfrg -frag AUTO_FORMAT.frg.fasta -pair AUTO_FORMAT.prs.fasta");
				}else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -ids $ids -numfrg $numfrg -frag AUTO_FORMAT.frg.fastq -pair AUTO_FORMAT.prs.fastq"); }
			
			}elsif( ($run_as_frag_bool == 1) && ($toggle == 1) ){
			
				my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
					
				my $numfrg; 
				if ($orig_fragments){ 
					$numfrg = `grep -c \'\>\' $orig_fragments`;
					chomp $numfrg;
				}else{
					$numfrg = 0;
				}
				
				if ($fastq_bool == 0){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -ids $ids -numfrg $numfrg -pair AUTO_FORMAT.frg.fasta");
				}else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -ids $ids -numfrg $numfrg -pair AUTO_FORMAT.frg.fastq -q");  }
					
				if ($fastq_bool == 0){ printl("Fragment only output provided as:\n\t+ fasta : NEATFREQ_OUT.all_fragments.fasta , NEATFREQ_OUT.all_mates.fasta\n"); }
				else{ printl("Fragment only output provided as:\n\t+ fasta :NEATFREQ_OUT.all_fragments.fasta , NEATFREQ_OUT.all_mates.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq , NEATFREQ_OUT.all_mates.fastq"); }
			
			}else{
				
				if ($fastq_bool == 0){ 
					if ( $input_status==2 ){
						# both prepare
						printl("Fragment only output provided as:\n\t+ fasta : $NF_prefix [..] X.fragments.fasta\n"); 
						printl("Paired-end only output provided as:\n\t+ fasta : $NF_prefix [..] X.pairs.fasta\n"); 
					}elsif ( $input_status==1 ){
						# pair prepare
						printl("Paired-end only output provided as:\n\t+ fasta : $NF_prefix [..] X.pairs.fasta\n"); 
					}elsif( $input_status==0 ){
						# fragments prepare
						printl("Fragment only output provided as:\n\t+ fasta : $NF_prefix [..] X.fragments.fasta\n"); 
					}
					
				}else{ 
				
					# UPDATE
					my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
					my $frg_ids = "$NF_prefix" . ".ALL.frags.ids.txt";
					my $prs_ids = "$NF_prefix" . ".ALL.pairs.ids.txt";
					if ( $input_status==2 ){
						# both exist
						
							# fragment convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $frg_ids -outfile NEATFREQ_OUT.all_fragments.fastq");
						printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
							# pair convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $prs_ids -outfile NEATFREQ_OUT.all_pairs.fastq");
						printl("Paired end output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.pairs.fastq\n"); 
						
					}elsif( $input_status==1 ){
						# pairs only
						
							# fragment convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $frg_ids -outfile NEATFREQ_OUT.all_fragments.fastq");
						printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
							# pair convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $prs_ids -outfile NEATFREQ_OUT.all_pairs.fastq");
						printl("Paired end output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.pairs.fastq\n"); 
						
					}elsif( $input_status==0 ){
						# fragments only
						
							# fragment convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $ids -outfile NEATFREQ_OUT.all_fragments.fastq");
						printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
						
					}
					
					#my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
					#runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile AUTO_FORMAT.frg.fastq -name $ids -outfile NEATFREQ_OUT.all_fragments.fastq");
					#printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
				
				}
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
			if ( $run_as_frag_bool == 0 && ($bin_extract_type eq "random") && $rpairs && (-s("$rpairs")) && $reads_file && (-s("$reads_file")) ){
					
				my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
				
				my $numfrg = `grep -c \'\>\' AUTO_FORMAT.frg.fasta`;
				chomp $numfrg;
				if ( (-s("AUTO_FORMAT.frg.fasta")) && (!($numfrg > 0)) ){
					printl("\n\nERROR : Failure to parse fasta output of NeatFreq.pl.  See NeatFreq_auto.log for details.\n\n");
					exit;
				}else{
					$numfrg = 0;
				}
				
				if ($fastq_bool == 0){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -N $NEW_NEATFREQ_INSTALL -ids $ids -numfrg $numfrg -frag AUTO_FORMAT.frg.fasta -pair AUTO_FORMAT.prs.fasta");
				}else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -N $NEW_NEATFREQ_INSTALL -ids $ids -numfrg $numfrg -frag AUTO_FORMAT.frg.fastq -pair AUTO_FORMAT.prs.fastq"); }
			
			}elsif( ($run_as_frag_bool == 1) && ($toggle == 1) ){
			
				my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
					
				my $numfrg; 
				if ($orig_fragments){ 
					$numfrg = `grep -c \'\>\' $orig_fragments`;
					chomp $numfrg;
				}else{
					$numfrg = 0;
				}
				
				if ($fastq_bool == 0){ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -N $NEW_NEATFREQ_INSTALL -ids $ids -numfrg $numfrg -pair AUTO_FORMAT.frg.fasta");
				}else{ runsys("perl $NEATFREQ_INSTALL/NeatFreq_postprocess_random.pl -N $NEW_NEATFREQ_INSTALL -ids $ids -numfrg $numfrg -pair AUTO_FORMAT.frg.fastq -q");  }
					
				if ($fastq_bool == 0){ printl("Fragment only output provided as:\n\t+ fasta : NEATFREQ_OUT.all_fragments.fasta , NEATFREQ_OUT.all_mates.fasta\n"); }
				else{ printl("Fragment only output provided as:\n\t+ fasta :NEATFREQ_OUT.all_fragments.fasta , NEATFREQ_OUT.all_mates.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq , NEATFREQ_OUT.all_mates.fastq"); }
			
			}else{
				
				if ($fastq_bool == 0){ 
					if ( $input_status==2 ){
						# both prepare
						printl("Fragment only output provided as:\n\t+ fasta : $NF_prefix [..] X.fragments.fasta\n"); 
						printl("Paired-end only output provided as:\n\t+ fasta : $NF_prefix [..] X.pairs.fasta\n"); 
					}elsif ( $input_status==1 ){
						# pair prepare
						printl("Paired-end only output provided as:\n\t+ fasta : $NF_prefix [..] X.pairs.fasta\n"); 
					}elsif( $input_status==0 ){
						# fragments prepare
						printl("Fragment only output provided as:\n\t+ fasta : $NF_prefix [..] X.fragments.fasta\n"); 
					}
						
				}else{ 
					
					# UPDATE
					my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
					my $frg_ids = "$NF_prefix" . ".ALL.frags.ids.txt";
					my $prs_ids = "$NF_prefix" . ".ALL.pairs.ids.txt";
					if ( $input_status==2 ){
						# both exist
						
							# fragment convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $frg_ids -outfile NEATFREQ_OUT.all_fragments.fastq");
						printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
							# pair convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $prs_ids -outfile NEATFREQ_OUT.all_pairs.fastq");
						printl("Paired end output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.pairs.fastq\n"); 
						
					}elsif( $input_status==1 ){
						# pairs only
						
							# fragment convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $frg_ids -outfile NEATFREQ_OUT.all_fragments.fastq");
						printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
							# pair convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $prs_ids -outfile NEATFREQ_OUT.all_pairs.fastq");
						printl("Paired end output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.pairs.fastq\n"); 
						
					}elsif( $input_status==0 ){
						# fragments only
						
							# fragment convert
						runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile allreads.fastq -name $ids -outfile NEATFREQ_OUT.all_fragments.fastq");
						printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
						
					}
					
					#my $ids = "$NF_prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
					#runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile AUTO_FORMAT.frg.fastq -name $ids -outfile NEATFREQ_OUT.all_fragments.fastq");
					#printl("Fragment only output provided as:\n\t+ fasta :  $NF_prefix [..] x.fasta\n\t+ fastq : NEATFREQ_OUT.all_fragments.fastq\n"); 
				
				}
			}
	}
		
	#CLEANUP
	printl("CLEANING UP!\n");
	system("rm -rf KMER_BINS_mp_frg KMER_BINS INITBINS_CDHITEST");
	my $outprefix = "$prefix" . ".NeatFreq_auto.LOG.txt";
	system ("mv ./run.LOG.txt $outprefix");
		
	close LOG;
}
