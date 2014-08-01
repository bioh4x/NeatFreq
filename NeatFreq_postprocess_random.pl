#!/usr/bin/perl -w
#INITIAL SETUP
# J. Craig Venter Institute
# NeatFreq ver 1.0.1
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)

use strict;
use File::Basename; #required for fileparse()\
use Getopt::Std; #better command line input method
use Getopt::Long; #better command line input method

#OPEN LOG
my $NEATFREQ_INSTALL = "/home/jmccorri/scripts/NF_prod_testing";
open LOG, "> ./NeatFreq_postprocess_random.LOG.txt";

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
	printl("NeatFreq_postprocess_random.pl\n");
	printl("\tInput Flags:\n");
	printl("\t-frag fragments.fast[a/1]\n\t-pair pairs.interleaved.fast[a/q]\n");
	printl("\t-q\ttoggles fasta or fastq mode (not auto-detected!)\n");
	printl("\t-ids ids.txt\tOne ID per line.\n\t-numfrg\tfragment count in original input (grep -c \'\>\' AUTO_FORMAT.frg.fasta)\n");
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
	my $id_file = shift;
	my $numfrg = shift;
	
	#FIX IDS
	runsys("cat $id_file | sort -n | uniq | awk \'{if (\$1 \<= $numfrg) print}\' > all_fragment_reads.ids.txt");
	runsys("cat $id_file | sort -n | uniq | awk \'{if (\$1 \> $numfrg) print}\' > all_mate_reads.ids.txt");


	#PAIRS - FIND MATES
	open (VALIDMATES, '> validmate.ids.txt');
	open (SINGLEMATES, '> matesasfragments.ids.txt');

	my $even_bool;
	if ($numfrg % 2 == 0){
		$even_bool = 1;
	}else{
		$even_bool = 0;
	}
	
	open(GC,"< all_mate_reads.ids.txt");
	my $id1 = -2; 
	my $id2 = -1;
	my $id1_hold;
	my $line;
	while (defined($line = <GC>)){ #get each line
		chomp $line;
		if ($line % 2 == $even_bool){ #toggled for even or odd offset!
			if (($id1 > $id2) && ($id2 != -1)){
				#front of mate found, does not match with a reverse id
				print SINGLEMATES "$id1\n";
			}
			#front of mate found first
			$id1 = "$line";
		}else{
			$id2 = "$line";
			if ( ($id1 == ($id2-1)) && ($id1 != -2) ){
				#back of mate found and makes second half of mate relationship
				#print 2-sided mate
				print VALIDMATES "$id1\n";
				print VALIDMATES "$id2\n";
			}else{
				#back of mate found, does not match up to previous forward id
				print SINGLEMATES "$id2\n";
			}
		}
	}
	if ( ($line != $id1) && ($id1 > $id2) && ($id2 != -1)){
		print SINGLEMATES "$id1\n";
	}
	close VALIDMATES;
	close SINGLEMATES;
	
	#FRAGMENTS - EXTRACT ALL
			
		# EXTRACT
		if ($fragfile && (-s($fragfile))){
			runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $fragfile -idlist all_fragment_reads.ids.txt -o NEATFREQ_OUT.fragsfromfrags");
		}
		runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $pairfile -idlist matesasfragments.ids.txt -o NEATFREQ_OUT.fragsfrommates");
		runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $pairfile -idlist validmate.ids.txt -o NEATFREQ_OUT.all_mates");
		
		#JOIN
		if (-s($fragfile)){
			runsys("cat NEATFREQ_OUT.fragsfromfrags_1.fasta > NEATFREQ_OUT.all_fragments.fasta");
		}
		runsys("cat NEATFREQ_OUT.fragsfrommates_1.fasta >> NEATFREQ_OUT.all_fragments.fasta");
		runsys("mv NEATFREQ_OUT.all_mates_1.fasta NEATFREQ_OUT.all_mates.fasta");
	
	#CLEAN UP
	runsys("rm all_fragment_reads.ids.txt");
	runsys("rm all_mate_reads.ids.txt");
	runsys("rm validmate.ids.txt");
	runsys("rm matesasfragments.ids.txt");
	runsys("rm NEATFREQ_OUT.fragsfromfrags_1.fasta");
	runsys("rm NEATFREQ_OUT.fragsfrommates_1.fasta");
	runsys("rm extractFasta.log");
}



########################################################################################################################
# fixfastq: 
########################################################################################################################
sub fixfastq {
	my $fragfile = shift;
	my $pairfile = shift;
	my $id_file = shift;
	my $numfrg = shift;
	
	#FIX IDS
	runsys("cat $id_file | sort -n | uniq | awk \'{if (\$1 \<= $numfrg) print}\' > all_fragment_reads.ids.txt");
	runsys("cat $id_file | sort -n | uniq | awk \'{if (\$1 \> $numfrg) print}\' > all_mate_reads.ids.txt");


	#PAIRS - FIND MATES
	open (VALIDMATES, '> validmate.ids.txt');
	open (SINGLEMATES, '> matesasfragments.ids.txt');

	my $even_bool;
	if ($numfrg % 2 == 0){
		$even_bool = 1;
	}else{
		$even_bool = 0;
	}
	
	open(GC,"< all_mate_reads.ids.txt");
	my $id1 = -2; 
	my $id2 = -1;
	my $id1_hold;
	my $line;
	while (defined($line = <GC>)){ #get each line
		chomp $line;
		if ($line % 2 == $even_bool){ #toggled for even or odd offset!
			if (($id1 > $id2) && ($id2 != -1)){
				#front of mate found, does not match with a reverse id
				print SINGLEMATES "$id1\n";
			}
			#front of mate found first
			$id1 = "$line";
		}else{
			$id2 = "$line";
			if ( ($id1 == ($id2-1)) && ($id1 != -2) ){
				#back of mate found and makes second half of mate relationship
				#print 2-sided mate
				print VALIDMATES "$id1\n";
				print VALIDMATES "$id2\n";
			}else{
				#back of mate found, does not match up to previous forward id
				print SINGLEMATES "$id2\n";
			}
		}
	}
	if ( ($line != $id1) && ($id1 > $id2) && ($id2 != -1)){
		print SINGLEMATES "$id1\n";
	}
	close VALIDMATES;
	close SINGLEMATES;
	
	#FRAGMENTS - EXTRACT ALL
			
		# EXTRACT
		#runsys("/usr/local/devel/BCIS/assembly/tools/extractFasta -i $fragfile -idlist all_fragment_reads.ids.txt -o NEATFREQ_OUT.fragsfromfrags");
		#runsys("/usr/local/devel/BCIS/assembly/tools/extractFasta -i $pairfile -idlist matesasfragments.ids.txt -o NEATFREQ_OUT.fragsfrommates");
		#runsys("/usr/local/devel/BCIS/assembly/tools/extractFasta -i $pairfile -idlist validmate.ids.txt -o NEATFREQ_OUT.all_mates");
		
		if ($fragfile && (-s($fragfile))){
			runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $fragfile -name all_fragment_reads.ids.txt -outfile NEATFREQ_OUT.fragsfromfrags_1.fastq");
		}
		runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $pairfile -name matesasfragments.ids.txt -outfile NEATFREQ_OUT.fragsfrommates_1.fastq");
		runsys("$NEATFREQ_INSTALL/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $pairfile -name validmate.ids.txt -outfile NEATFREQ_OUT.all_mates.fastq");
		
		#JOIN
		runsys("cat NEATFREQ_OUT.fragsfromfrags_1.fastq > NEATFREQ_OUT.all_fragments.fastq");
		runsys("cat NEATFREQ_OUT.fragsfrommates_1.fastq >> NEATFREQ_OUT.all_fragments.fastq");
		#runsys("mv NEATFREQ_OUT.all_mates_1.fastq NEATFREQ_OUT.all_mates.fasta");
	
	#CLEAN UP
	runsys("rm all_fragment_reads.ids.txt");
	runsys("rm all_mate_reads.ids.txt");
	#runsys("rm validmate.ids.txt");
	runsys("rm matesasfragments.ids.txt");
	runsys("rm NEATFREQ_OUT.fragsfromfrags_1.fastq");
	runsys("rm NEATFREQ_OUT.fragsfrommates_1.fastq");
}



########################################################################################################################
#       MAIN                                                                                                           #
########################################################################################################################
MAIN : {
#command line input vars
	my %Opts;
	my ($id_file, $numfrg, $fragfile, $pairfile, $fastq_status, $NEW_NEATFREQ_INSTALL);
	
#pull in command line input
	my $status = GetOptions(\%Opts, "help", "h", "q", 'ids=s'=> \$id_file, 'numfrg=s'=> \$numfrg, 'N=s'=> \$NEW_NEATFREQ_INSTALL, 'frag=s'=> \$fragfile, 'pair=s'=> \$pairfile);
	if ( exists $Opts{help} || exists $Opts{h}){ outhelp(); }
	if ( exists $Opts{q} ){ 
		printl("FASTQ MODE enabled\n");
		$fastq_status = 1; }
	else{ $fastq_status = 0; }
	if ($numfrg && $numfrg > 0){
		printl("\nUSING INPUT FRAMENT COUNT = $numfrg\n");
	}elsif($numfrg == 0){
		printl("\nUSING NULL FRAGMENT COUNT = $numfrg\n");
		
	}else{
		printl("\nCAN'T FIND = $numfrg\n");
		outhelp();
	}
	if (($NEW_NEATFREQ_INSTALL) && ($NEW_NEATFREQ_INSTALL ne $NEATFREQ_INSTALL)){
		printl("Forcing use of user-selected NeatFreq install directory : $NEW_NEATFREQ_INSTALL\n");
		if (!(-s("$NEW_NEATFREQ_INSTALL/NeatFreq.pl"))){
			printl("\n\nERROR! : No NeatFreq install found in suggested install directory.  \n\nExiting...\n");
			exit;
		}
		$NEATFREQ_INSTALL = $NEW_NEATFREQ_INSTALL;
	}
	if ( $id_file ){ 
		if (-s ("$id_file")) { 
			printl("\nUSING INPUT ID FILE : $id_file\n"); 
		}
		else{ 
			printl("\nPROBLEM FOUND : File $id_file does not exist or has 0 byte size.\nExiting...");
			exit;
		}  
	}
	if ( $fragfile ){ 
		if (-s ("$fragfile")) { 
			printl("\nUSING INPUT CAS FILE : $fragfile\n"); 
		}
		else{ 
			#printl("\nPROBLEM FOUND : File $fragfile does not exist or has 0 byte size.\nExiting...");
			printl("\nWARNING : File $fragfile does not exist or has 0 byte size.\n");
			#exit;
		}  
	}
	if ( $pairfile ){ 
		if (-s ("$pairfile")) { 
			printl("\nUSING INPUT CAS FILE : $pairfile\n"); 
		}
		else{ 
			printl("\nPROBLEM FOUND : File $pairfile does not exist or has 0 byte size.\nExiting...");
			exit;
		}  
	}
	
#Fix 'em
	if ($fastq_status == 0){
		fixfasta($fragfile, $pairfile, $id_file, $numfrg);
	}else{
		fixfastq($fragfile, $pairfile, $id_file, $numfrg);
	}
	
	exit;
}
