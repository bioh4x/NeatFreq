#!/usr/bin/perl -w
use strict;
use File::Basename; #required for fileparse()
use Getopt::Std; #better command line input method
use Getopt::Long; #better command line input method

#INSTALL LOCATIONS (1st party)
my $NEATFREQ_INSTALL = "/usr/local/devel/BCIS/assembly/tools/NeatFreq";
my $AUTOMATON_INSTALL = "$NEATFREQ_INSTALL" . "/assembly_automaton/lib";
#INSTALL LOCATIONS (3rd party)
my $MIRA_INSTALL = "$NEATFREQ_INSTALL" . "/lib";
my $CLC_NGS_INSTALL = "/usr/local/packages/clc-ngs-cell";
my $CUTADAPT_INSTALL = "/usr/local/devel/BCIS/external_software/cutadapt-1.0";
my $QUAKE_INSTALL = "/usr/local/devel/BCIS/external_software/quake-0.3.0";
my $PYTHON_INSTALL = "/usr/local/bin"; # must be configured for Allpaths Error Correction
my $PERL_INSTALL = "/usr/local/bin";
my $FASTQ_UTIL_INSTALL = "/usr/local/devel/BCIS/pratap/bin/utilities/brentp-bio-playground-38c7970/reads-utils";
my $DUST_INSTALL = "/local/platform/bin";
my $RNNOTATOR_INSTALL = "/usr/local/packages/Rnnotator";
my $CELERA_WGS_INSTALL = "/usr/local/packages/wgs/Linux-amd64";
my $VELVET_INSTALL = "/usr/local/packages/velvet";
my $SEQ454_INSTALL = "/usr/local/packages/seq454"; # also known as Newbler runAssembly install

#GLOBALS
my %Opts;
my $pwd = `pwd`; chomp $pwd;
my $DIR = $pwd;
my @adapters;
my $totalDUSTchanges = 0;
my $no_assemb_switch;

#OPEN LOG
if (-e("./assembly_automaton.LOG.txt")){ system("rm ./assembly_automaton.LOG.txt"); }
system("touch ./assembly_automaton.LOG.txt");
open LOG, ">> ./assembly_automaton.LOG.txt";

################################################################################################################################################################################################################################################
# printl : 
# print anything once to STDOUT and once to the LOG file -- pretty useful, eh?
################################################################################################################################################################################################################################################
sub printl {
	my $pr = shift;
	#run STDOUT print
	print $pr;
	#run LOG print (when appropriate)
	print LOG $pr;
}

################################################################################################################################################################################################################################################
# runsys : 
################################################################################################################################################################################################################################################
sub runsys {
	my $cmd = shift;
	printl("++ RUN CMD : $cmd\n");
	system("$cmd");
}

################################################################################################################################################################################################################################################
# outhelp : 
# print anything once to STDOUT and once to the LOG file -- pretty useful, eh?
################################################################################################################################################################################################################################################
sub outhelp {
	print "\n          _____             \n";
	print "         /_____\\           \n";
	print "    ____[\\`---'/]____      \n";
	print "   /\\ #\\ \\_____/ /# /\\     ASSEMBLY AUTOMATON\n";
	print "  /  \\# \\_.---._/ #/  \\    ver 1.0 - 2/2012\n";
	print " /   /|\\  |   |  /|\\   \\   jmccorri\@jcvi.org\n";
	print "/___/ | | |   | | | \\___\\  \n";
	print "|  |  | | |---| | |  |  |  Takes in sff,fasta,fastq representations\n";
	print "|__|  \\_| |_#_| |_/  |__|  of 454 and illumina sequencing and runs\n";
	print "//\\\\  <\\ _//^\\\\_ />  //\\\\  pre-processing and assembly pipelines.\n";
	print "\\||/  |\\//// \\\\\\\\/|  \\||/  \n";
	print "      |   |   |   |        \n";
	print "      |___|   |___|        \n";
	print "      /   \\   /   \\        \n";
	print "     |_____| |_____|       \n\n";
	print "EXAMPLE USAGE : \nassembly_automaton.pl -r -u -c -z -o ./OUT_DIR -g auto -p PREFIX -SFFfragment in1.sff in2.sff -ILLUMpair interleaved.1and2.fastq -ILLUMinsert 260\n\n";
	print "Required input:\n";
	print "\t-o <Output Dir>\n";
	print "\t-g <Genome Size>\n";
	print "Sequence Data input:\n";
	print "\t-SFFfragment <input1.sff -SFFfragment input2.sff ...>\n";
	print "\t-SFFpair <input1.sff -SFFpair input2.sff ...>\n";
	print "\t-SFFinsert <sff insert size>\t\tStd. dev. auto set to 10% insert size.\n";
	print "\t-ILLUMfragment <input1.fastq -ILLUMfragment input2.fastq ...>\n";
	print "\t-ILLUMpair <input1.fastq -ILLUMpair input2.fastq ...>\n";
	print "\t\t(Must be interleaved.)\n";
	print "\t\t($AUTOMATON_INSTALL/shuffleSequences_fastq.pl)\n";
	print "\t-ILLUMinsert <illum. insert size>\tStd. dev. auto set to 10% insert size.\n";
	print "Optional input:\n";
	print "\t-r <quake/allpaths>: Turn on Read Correction [illumina-only] and specify program you'd like to use. (for high cov)\n";
	print "\t-u : Turn on Uniqueness/De-duplication Filtering. [illumina-only] (for high cov, does not work with -r allpaths)\n";
	print "\t-d : Turn on DUST low complexity masking.\n";
	#print "\t-c : Turn on Random Coverage Filtering Stage. (for high unique cov)\n";
	#print "\t-z : Turn on Post-Assembly Chaff Analysis and Filtering. (for uneven cov)\n";
	print "\t-adapter <adapter sequence> : Cut adapters from sequence.  Can be used multiple times\n";
	print "\t-p <Output Prefix>\t\tDefault = OUT\n";
	print "\t-f <Illumina QV Offset>\t\tDefault = 0\n";
	print "\t-precontam <Contaminant Fasta>\t Default = <none>\n";
	print "\t-contam <Contaminant Fasta>\t Default = <none>\n";
	print "\t\t(Merge multiple contaminants into single fasta file.)\n";
	print "\t-help or -h :\tShow this help output.\n\n";
	print "Filter pipeline order:\n";
	print "Pre-sub contam Check -> Illumina Read Correction -> [Uniqueness Check] -> Low Complexity Mask -> Illum. QV trim -> [*]\n";
	print "[*] -> Contam Check -> Adapter Removal -> Assembly -> Post-Assembly Processing -> Clean Suspected Chaff\n\nExiting.\n";
	
	close LOG;
	exit;
}


################################################################################################################################################################################################################################################
# reverse_complement
################################################################################################################################################################################################################################################
sub reverse_complement {
	my $dna = shift;
	
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	
	print "Reverse complementing [ $dna ] -> [ $revcomp ]\n";
	
	return $revcomp;
}

################################################################################################################################################################################################################################################
# sequence_info
################################################################################################################################################################################################################################################
sub sequence_info {
	my @inputs = @_;
	my @return;
	
	my $c = 0;
	my @readsizeX;
	my @readcountX;
	foreach my $fasta (@inputs){
		runsys("$CLC_NGS_INSTALL/sequence_info $fasta > $fasta.sequence_info");
		
		my $tot = `grep Number $fasta.sequence_info`;
		chomp $tot;
		my @xs = split(/\s+/, $tot);
		$tot = $xs[3];
		$readcountX[$c] = $xs[3];	
		#printl("NUMBER OF READS : $tot\n");
		$return[0] += $tot;
		
		$tot = `grep Total $fasta.sequence_info`;
		chomp $tot;
		@xs = split(/\s+/, $tot);
		$tot = $xs[2];
		#printl("BASES IN READS : $tot\n");
		$return[1] += $tot;
		
		my $readsize = `grep Average $fasta.sequence_info`;
		chomp $readsize;
		@xs = split(/\s+/, $readsize);
		$readsizeX[$c] = $xs[2];	
	}
	
	my $tmpavg=0;
	my $d=0;
	my $runtot = 0;
	foreach(@readsizeX){
		if ($_){
			$tmpavg = $_ * $readcountX[$d];
			$runtot += $readcountX[$d];
			$d++;
		}
	}
	if ($runtot){
		$return[2] = $tmpavg/$runtot;
	}
	
	return @return;
}


################################################################################################################################################################################################################################################
# sequence_info
################################################################################################################################################################################################################################################
sub summarize_reads{
	#my %dataset_hash = %{$_[0]};
	my %dataset_hash = @_;
	my @tmp1; my @tmp2; my @tmp3; my @tmp4; my @tmp5;
	my $c1=0; my $c2=0; my $c3=0; my $c4=0; my $c5=0;
	if ( @{$dataset_hash{0}} ){ 
		@tmp1 = sequence_info(@{$dataset_hash{0}}); 
		$c1 = @{$dataset_hash{0}};
	}
	if ( @{$dataset_hash{1}} ){ 
		@tmp2 = sequence_info(@{$dataset_hash{1}}); 
		$c2 = @{$dataset_hash{1}};
	}
	if ( @{$dataset_hash{2}} ){ 
		printl("\n\nEVALUATING ILLUMINA FRAGMENTS:\n");
		@tmp3 = sequence_info(@{$dataset_hash{2}});
		$c3 = @{$dataset_hash{2}}; 
		printl("END EVALUATION OF ILLUMINA FRAGMENTS\n\n\n");
	}
	if ( @{$dataset_hash{3}} ){ 
		@tmp4 = sequence_info(@{$dataset_hash{3}}); 
		$c4 = @{$dataset_hash{3}};
	}
	if ( @{$dataset_hash{4}} ){ 
		@tmp5 = sequence_info(@{$dataset_hash{4}}); 
		$c5 = @{$dataset_hash{4}};
	}
	if (!($tmp1[1])){ $tmp1[0] = 0; $tmp1[1] = 0; }
	if (!($tmp2[1])){ $tmp2[0] = 0; $tmp2[1] = 0; }
	if (!($tmp3[1])){ $tmp3[0] = 0; $tmp3[1] = 0; }
	if (!($tmp4[1])){ $tmp4[0] = 0; $tmp4[1] = 0; }
	if (!($tmp5[1])){ $tmp5[0] = 0; $tmp5[1] = 0; }
	
	my $in_tot = $c1 + $c2 + $c3 + $c4 + $c5;
	my $read_count = $tmp1[0] + $tmp2[0] + $tmp3[0] + $tmp4[0] + $tmp5[0];
	my $bp_count = $tmp1[1] + $tmp2[1] + $tmp3[1] + $tmp4[1] + $tmp5[1];
	my $tmp_avg; my $tmp;
	printl("\t\t\t\t\tFiles\tPerc.\tNum.Reads\tNum.BPs\t\tPerc.BPs\tAvg.Read.Len\n");
	if ( ($tmp1[1]) && ($tmp1[1] > 0) ){
		$tmp_avg = sprintf( "%.2f", (($tmp1[0]/$read_count)*100) ) ;
		$tmp = sprintf( "%.2f", (($c1/$in_tot)*100) );
		printl("Total fragment-only 454 :\t\t$c1\t$tmp%\t$tmp1[0] reads\t$tmp1[1] bp\t$tmp_avg%\t\t$tmp1[2] bp\n");
	}else{
		printl("Total fragment-only 454 :\t\t0\t0.00 %\n");
	}
	if ( ($tmp2[1]) && ($tmp2[1] > 0) ){
		$tmp_avg = sprintf( "%.2f", (($tmp2[0]/$read_count)*100) ) ;
		$tmp = sprintf( "%.2f", (($c2/$in_tot)*100) );
		printl("Total paired-end 454 :\t\t\t$c2\t$tmp%\t$tmp2[0] reads\t$tmp2[1] bp\t$tmp_avg%\t\t$tmp2[2] bp\n");
	}else{
		printl("Total paired-end 454 :\t\t\t0\t0.00 %\n");
	}
	if ( ($tmp3[1]) && ($tmp3[1] > 0) ){
		$tmp_avg = sprintf( "%.2f", (($tmp3[0]/$read_count)*100) ) ;
		$tmp = sprintf( "%.2f", (($c3/$in_tot)*100) );
		printl("Total fragment-only Illum :\t\t$c3\t$tmp%\t$tmp3[0] reads\t$tmp3[1] bp\t$tmp_avg%\t\t$tmp3[2] bp\n");
	}else{
		printl("Total fragment-only Illum :\t\t0\t0.00 %\n");
	}
	if ( ($tmp4[1]) && ($tmp4[1] > 0) ){
		$tmp_avg = sprintf( "%.2f", (($tmp4[0]/$read_count)*100) ) ;
		$tmp = sprintf( "%.2f", (($c4/$in_tot)*100) );
		printl("Total paired-end Illum :\t\t$c4\t$tmp%\t$tmp4[0] reads\t$tmp4[1] bp\t$tmp_avg%\t\t$tmp4[2] bp\n");
	}else{
		printl("Total paired-end Illum :\t\t0\t0.00 %\n");
	}
	if ( ($tmp5[1]) && ($tmp5[1] > 0) ){
		$tmp_avg = sprintf( "%.2f", (($tmp5[0]/$read_count)*100) ) ;
		$tmp = sprintf( "%.2f", (($c5/$in_tot)*100) );
		printl("Total fragment-only fasta :\t\t$c5\t$tmp%\t$tmp5[0] reads\t$tmp5[1] bp\t$tmp_avg%\t\t$tmp5[2] bp\n");
	}else{
		#printl("Total fragment-only fasta :\t\t0\t0.00 %\n");
	}
	printl("TOTAL:\t\t\t\t\t$in_tot\t-\t$read_count reads\t$bp_count bp\n\n");
	
	return $bp_count;
}


################################################################################################################################################################################################################################################
# contam_rm
################################################################################################################################################################################################################################################
sub contam_rm {
	my $contam_ref = shift;
	my $is_presub = shift; #1=presub #0=contam
	my $format = shift; #sff fasta fastq
	my $is_pairs = shift;
	my @query_seqs = @_;
	my @contam_rm;
	my @contam_rm_newfragments;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	if ($is_presub == 1){ 
		runsys("mkdir $DIR/presub_contam"); 
		chdir("$DIR/presub_contam") or die "Could not change directory to $DIR/presub_contam\n"; 
	}else{ 
		runsys("mkdir $DIR/contam"); 
		chdir("$DIR/contam") or die "Could not change directory to $DIR/contam\n"; 
	}
	foreach my $seq (@query_seqs){
		my($filename, $directory, $suffix) = fileparse("$seq");
		
			#CORRECT QUAL ID LINE?
			my $tmp_in = `head -n 3 $seq | tail -n 1`; chomp $tmp_in;
			my $tmp_out = "$DIR" . "/" . "$filename" . ".IDfix.fastq";
			if ($tmp_in eq "+"){ 
				#printl("Fixing ID headers on quality lines of fastq file...\n");
				runsys("$AUTOMATON_INSTALL/fastq_ID_fixer.pl -i $seq > $tmp_out");			
				if (-s($tmp_out)){
					$seq = $tmp_out;
					($filename, $directory, $suffix) = fileparse("$seq");
				}else{
					printl("\nERROR during initial QV id line fix in contamination removal..\n");
					exit;
				}
			}
		
		my $tmpfile = "$filename" . ".fasta";
		my $tmpcas = "$filename" . ".vs_contam.cas";
		my $tmptable = "$filename" . ".vs_contam.assembly_table.txt";
		my $tmpcontamids = "$filename" . ".contam_ids.txt";
		#align to ref
		if ($format eq 'sff'){
			runsys("$SEQ454_INSTALL/bin/sffinfo -s $seq > $tmpfile");
			runsys("ln -s $CLC_NGS_INSTALL/license.properties .");
			runsys("$CLC_NGS_INSTALL/clc_ref_assemble_long -o $tmpcas -l 0.4 -s 0.95 -q $tmpfile -d $contam_ref");
		}else{
			runsys("ln -s $CLC_NGS_INSTALL/license.properties .");
			runsys("$CLC_NGS_INSTALL/clc_ref_assemble_long -o $tmpcas -l 0.4 -s 0.95 -q $seq -d $contam_ref");
		}
		#build and parse assembly table
		runsys("$CLC_NGS_INSTALL/assembly_table -n $tmpcas > $tmptable");
		#
		#extract
		if ($format eq 'sff'){
			runsys("cat $tmptable | awk \'{if (\$10 > -1) print \$2}\' | sort > $tmpcontamids");
			printl("Extracting contaminated read ID's... ( $filename )\n");
			my $new_sff = "$filename" . ".contamRM.sff";
			runsys("$SEQ454_INSTALL/bin/sfffile -o $new_sff -e $tmpcontamids $seq");
			if (-s($new_sff)){
				my $tmp_pwd = `pwd`; chomp $tmp_pwd;
				$new_sff = "$tmp_pwd" . "/" . "$new_sff";
				push(@contam_rm, $new_sff);
			}else{
				printl("\nERROR during contaminant removal on $filename using $contam_ref -- See assembly_automaton.log for information on last process.\n");
				exit;
			}
		}else{
			runsys("cat $tmptable | awk \'{if (\$10 < 0) print \$2}\' | sort > $tmpcontamids");
			#EXTRACT FRAGMENTS VS MATES
			my $outfile;
			if ($is_pairs == 1){
				if (-s("$tmpcontamids")){
					runsys("$AUTOMATON_INSTALL/extract_fastq_as_mates_frags.pl $tmpcontamids OUT")
				}else{ printl("+ ERROR : No file $tmpcontamids found during contamination removal. +\nExiting..."); exit; }
				if ((-s("./OUT.mate_ids.txt"))){ 
					if ($is_presub == 1){ $outfile = "$DIR/presub_contam/" . "$filename" . ".contamRM_mates.fastq"; }
					else{ $outfile = "$DIR/contam/" . "$filename" . ".contamRM_mates.fastq"}
					runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $seq -name ./OUT.mate_ids.txt -outfile $outfile");
					if (-s("$outfile")){ push(@contam_rm, "$outfile"); }
					else{ printl("+ ERROR : Fastq extraction failed to build file [ $outfile ] +\nExiting...\n"); exit; }
				}else{ 
					printl("+ Adapter removal file [ ./OUT.mate_ids.txt ] not found. See logs. +\nExiting...\n"); 
					exit; 
				}
				if ((-s("./OUT.fragment_ids.txt"))){ 
					if ($is_presub == 1){ $outfile = "$DIR/presub_contam/" . "$filename" . ".contamRM_fragments.fastq"; }
					else{ $outfile = "$DIR/contam/" . "$filename" . ".contamRM_fragments.fastq"}
					runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $seq -name ./OUT.fragment_ids.txt -outfile $outfile");
					if (-s("$outfile")){ push(@contam_rm_newfragments, "$outfile"); }
					else{ printl("+ ERROR : Fastq extraction failed to build file [ $outfile ] +\nExiting...\n"); exit; }
				}else{ 
					printl("+ Adapter removal file [ ./OUT.fragment_ids.txt ] not found. See logs. +\n"); 
				}
			}else{
				if (-s("$tmpcontamids")){
					if ($is_presub == 1){ $outfile = "$DIR/presub_contam/" . "$filename" . ".contamRM_fragONLY.fastq"; }
					else{ $outfile = "$DIR/contam/" . "$filename" . ".contamRM_fragONLY.fastq"}
					runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $seq -name $tmpcontamids -outfile $outfile");
					if (-s("$outfile")){
						push(@contam_rm, "$outfile");
					}else{ printl("+ ERROR : Fastq extraction failed to build file [ $outfile ] +\nExiting...\n"); exit; }
				}else{ printl("+ ERROR : No file $tmpcontamids found during contamination removal. +\nExiting..."); exit; }
			}
			
			
		}
	}
	
	#SUMMARY
	printl("\n\nContaminant removal run complete.\n-> Output files parsed from pairs:\n");
	my $z = 1;
	foreach(@contam_rm){
		print ("$z - $_\n");
		$z++;
	}
	printl("-> From fragments only:\n");
	foreach (@contam_rm_newfragments){
		printl("$z - $_\n");
		$z++;
	}
	print "\n\n";
	
	my %new_hash;
	$new_hash{1} = [@contam_rm];
	$new_hash{2} = [@contam_rm_newfragments];
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return %new_hash;

}


################################################################################################################################################################################################################################################
# adapter_rm
################################################################################################################################################################################################################################################
sub adapter_rm {
	my $SFFinsert = shift;
	my $type = shift;
	my $is_pairs = shift; #0 or 1
	my $seq_count = shift;
	my @sequences = @_;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/adapter_rm"); 
	chdir("$DIR/adapter_rm") or die "Could not change directory to $DIR/presub_contam\n"; 
	
	my @new_sequences = ();
	my @newer_sequences_newfragments =();
	my $mock_count = 0;
	
	#CONVERT TO ILLUMINA IF SFF
	foreach my $sequence (@sequences){
		if ($type eq "sff"){ 
			if ($is_pairs == 0){
				my $libname = "frag" . "$mock_count";
				my $mock_out = "$libname" . ".mock.frg";
				my $convert_cmd = "$CELERA_WGS_INSTALL/bin/sffToCA";
				#$convert_cmd .= " -trim chop -libraryname $libname -output $mock_out $_";
				$convert_cmd .= " -trim chop -libraryname $libname -output $mock_out $sequence";
				runsys("$convert_cmd");
				if (-s("$mock_out")){
					my $libstore = "$libname" . "Store";
					runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
					if (-e("$libstore")){

						my $new_out = "$DIR/adapter_rm/" . "$libname" . ".unmated.fastq";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfastq $libname $libstore");
						if (-s($new_out)){
							push(@new_sequences, "$new_out");
						}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
				}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
				$mock_count++
			}else{
				my $libname = "frag" . "$mock_count";
				my $mock_out = "$libname" . ".mock.frg";
				my $convert_cmd = "$CELERA_WGS_INSTALL/bin/sffToCA";
				my $SFFsd = sprintf( "%.0f", ($SFFinsert/10));
				$convert_cmd .= " -linker titanium -insertsize $SFFinsert $SFFsd";
				#$convert_cmd .= " -trim chop -libraryname $libname -output $mock_out $_";
				$convert_cmd .= " -trim chop -libraryname $libname -output $mock_out $sequence";
				runsys("$convert_cmd");
				if (-s("$mock_out")){
					my $libstore = "$libname" . "Store";
					runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
					if (-e("$libstore")){

						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfastq $libname $libstore");
						my $new_out = "$DIR/adapter_rm/" . "$libname" . ".paired.fastq";
						if (-s($new_out)){
							push(@new_sequences, "$new_out");
						}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						$new_out = "$DIR/adapter_rm/" . "$libname" . ".unmated.fastq";
						if (-s($new_out)){
							push(@new_sequences, "$new_out");
						}else{ printl("Warning : no fragments output after adapter trimming on paired end input file.\n"); exit; }
					}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
				}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
				$mock_count++;
			}
		}else{ push(@new_sequences, "$sequence"); }
	}
	
	
	#NOW RUN CHECK ON ALL ILLUMINA DATA!
	my @newer_sequences;
	my @newer_sequences_newfragments;
	foreach my $seq (@new_sequences){
		my($filename, $directory, $suffix) = fileparse("$seq");
		foreach my $ad (@adapters){
			my $rev_ad = reverse_complement($ad);
			my $front_check_ad = "$ad" . "NNNNNN";
			my $back_check_ad = "NNNNNNN" . "$rev_ad";
			
			#cuts
			my $front_out = "$DIR/adapter_rm/FRONT_CUT." . "$filename";
			my $report = "$DIR/adapter_rm/FRONT_CUT." . "$filename" . "cutAdaptLog.txt";
			runsys("python $CUTADAPT_INSTALL/cutadapt -g $front_check_ad $seq -o $front_out > $report");
				if (!(-s("$front_out"))){ printl("+ Adapter removal file [ $front_out ] not found. See logs. +\nExiting...\n"); exit; }
			my $back_out = "$DIR/adapter_rm/BACK_CUT." . "$filename";
			$report = "$DIR/adapter_rm/BACK_CUT." . "$filename" . "cutAdaptLog.txt";
			runsys("python $CUTADAPT_INSTALL/cutadapt -a $back_check_ad $front_out -o $back_out > $report");
				if (!(-s("$back_out"))){ printl("+ Adapter removal file [ $back_out ] not found. See logs. +\nExiting...\n"); exit; }
			my $mid_out = "$DIR/adapter_rm/MID_CUT." . "$filename";
			$report = "$DIR/adapter_rm/MID_CUT." . "$filename" . "cutAdaptLog.txt";
			runsys("python $CUTADAPT_INSTALL/cutadapt -o 20 -b $ad -b $rev_ad --discard-trimmed $back_out -o $mid_out > $report");
				if (!(-s("$mid_out"))){ printl("+ Adapter removal file [ $mid_out ] not found. See logs. +\nExiting...\n"); exit; }
			if ($is_pairs == 1){
				runsys("cat $mid_out | grep \'\@\' | tr \'\@\' \' \' | awk \'{print \$1}\' > ids.txt");
					if (!(-s("./ids.txt"))){ printl("+ Adapter removal file [ ids.txt ] not found. See logs. +\nExiting...\n"); exit; }
				runsys("$AUTOMATON_INSTALL/extract_fastq_as_mates_frags.pl ids.txt OUT");
					if ((-s("./OUT.mate_ids.txt"))){ 
						my $outfile = "$seq" . ".contamRM_mates.fastq";
						runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $mid_out -name ./OUT.mate_ids.txt -outfile $outfile");
						if (-s("$outfile")){ push(@newer_sequences, "$outfile"); }
						else{ printl("+ ERROR : Fastq extraction failed to build file [ $outfile ] +\nExiting...\n"); exit; }
					}else{ 
						printl("+ Adapter removal file [ ./OUT.mate_ids.txt ] not found. See logs. +\nExiting...\n"); 
						exit; 
					}
					if ((-s("./OUT.fragment_ids.txt"))){ 
						my $outfile = "$DIR/adapter_rm/" . "$filename" . ".contamRM_fragments.fastq";
						runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $mid_out -name ./OUT.fragment_ids.txt -outfile $outfile");
						if (-s("$outfile")){ push(@newer_sequences_newfragments, "$outfile"); }
						else{ printl("+ ERROR : Fastq extraction failed to build file [ $outfile ] +\nExiting...\n"); exit; }
					}else{ printl("+ Adapter removal file [ ./OUT.fragment_ids.txt ] not found. See logs. +\n"); 
						exit; 
					}
			}else{
				push(@newer_sequences, "$mid_out");
			}
		}
	}
	
	# RETURN -- UPDATE MAIN TO ACCEPT ARRAY!!!!
	my %new_hash;
	$new_hash{1} = [@newer_sequences];
	$new_hash{2} = [@newer_sequences_newfragments];
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return %new_hash;
	#return @newer_sequences;
	
}


################################################################################################################################################################################################################################################
# read_correct_quake
################################################################################################################################################################################################################################################
sub read_correct_quake {
	my $offset = shift;
	my @query_seqs = @_;
	my @corrected_seqs;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/read_correct"); 
	chdir("$DIR/read_correct") or die "Could not change directory to $DIR/read_correct\n"; 
	
	if (-e("./fastq_files.txt")){
		runsys("rm ./fastq_files.txt");
	}
	runsys("touch ./fastq_files.txt");
	foreach my $seq (@query_seqs){
		runsys("echo $seq >> ./fastq_files.txt")
	}
	
	runsys("$QUAKE_INSTALL/bin/quake.py --no_jelly -f fastq_files.txt -k 19 -q $offset --log > quakepy.log");
	if (-s("fastq_files.txt.qcts")){
		runsys("$QUAKE_INSTALL/bin/correct -f fastq_files.txt -k 19 -m fastq_files.txt.qcts -c 12 -q $offset --log");
	}
	else{
		printl("\n\nQUAKE RUN FAILURE! NO QCTS FILE FOUND.\nExiting...\n");
	}
	
	#IDENTIFY OUTPUT FILE
	printl("Identifying corrected output...\n");
	if (-s ("./corrected_reads.txt")){ runsys("rm ./corrected_reads.txt"); }
	foreach my $seq (@query_seqs){
		my($filename, $directory, $suffix) = fileparse("$seq");
		my @zz = split(/./, $filename);
		my $search = "$directory";
		foreach my $z (@zz){
			if ($z ne "$suffix"){ $search .= "$z"; }
		}
		$search .= "*cor*";
		runsys("touch corrected_reads.txt");
		my $tmp = `ls $search`; chomp $tmp;
		if ($tmp){ 
			if (-s("$tmp")){
				my $newtmp =  "$DIR/read_correct/" . "$filename";
				runsys("mv $tmp $newtmp");
				push(@corrected_seqs, "$newtmp");
			}else{
				printl("+ File not found with expected name and location : $tmp\nSee logs for details.\nExiting...\n");
				exit;
			}
		}else{
			printl("+ No corrected read output found for library $filename. See logs.\nExiting...\n");
			exit;
		}
		
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return @corrected_seqs;
}


################################################################################################################################################################################################################################################
# read_correct_allpaths
################################################################################################################################################################################################################################################
sub read_correct_allpaths {
	my $offset = shift;
	my $ILLUMinsert = shift;
	my $FRAGMENTONLY = shift;
	my @query_seqs = @_;
	my @corrected_seqs;

	#DEBUG
	if ($FRAGMENTONLY == 1){
		$ILLUMinsert = 270;
	}
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/read_correct"); 
	chdir("$DIR/read_correct") or die "Could not change directory to $DIR/read_correct\n"; 
	
	my $correct_me; my $phred64;
	LINE: foreach my $seq (@query_seqs){
	#foreach my $seq (@query_seqs){
		
		#if cshell
		my($filename, $directory, $suffix) = fileparse("$seq");
		
		#qv check and move files
		if ($suffix ne "fastq"){
			$correct_me = "./" . "PE-" . "$ILLUMinsert" . "." . "$filename" . ".fastq";
		}else{
			$correct_me = "./" . "PE-" . "$ILLUMinsert" . "." . "$filename";
		}
		if ($offset == 64){
			$phred64 = 1;
			runsys("ln -s $seq $correct_me")
		}elsif ($offset == 33){
			$phred64 = 0;
			runsys("ln -s $seq $correct_me")
		}else{
			#modify qv
			#?
		}
		
		#correct
			runsys("$PYTHON_INSTALL/python -m jcvi.assembly.preprocess correct $correct_me");
		
		#find and move output
		my $expected_pair_1 = "$DIR/read_correct/" . "PE-" . "$ILLUMinsert" . ".1.corr.fastq";
		my $expected_pair_2 = "$DIR/read_correct/" . "PE-" . "$ILLUMinsert" . ".2.corr.fastq";
		
		
		#FAIL CASE
		if ( !(-s("$DIR/read_correct/frag_read_corr.corr.fastq")) && !(-s("$expected_pair_1")) && !(-s("expected_pair_2")) ){
			printl("+ ERROR : NO ALLPATHS OUTPUT FOUND IN $DIR/read_correct: +\n./frag_read_corr.corr.fastq\n$expected_pair_1\n$expected_pair_2\n\nExiting...");
			
			my $inputcount = `$CLC_NGS_INSTALL/sequence_info $correct_me | grep \'Number\' | awk \'{print \$4}\'`;
			chomp $inputcount;
			if ($inputcount < 2500){
				printl("+ WARNING : ALLPATHS parsing failed for small input file ( $inputcount reads) $filename.  Passing as corrected. +\n");
				push(@corrected_seqs, "$seq");
					#CLEANUP DIR
					printl("Cleaning read_correct directory...\n");
					runsys("rm -rf ./data");
					runsys("rm ./in_*");
				next; #move to next sequence
			}else{
				printl("+ ERROR : ALLPATHS parsing failed for large input file ( $inputcount reads) $filename.  Exiting. +\n");
				exit;	
			}
		}
			
		#FRAGMENT OUTPUT
		if (-s("./frag_reads_corr.corr.fastq")){
			#file found and contains sequence
			my $newname = "$DIR/read_correct/" . "$filename" . ".read_corr.0.fastq";
			runsys("mv ./frag_reads_corr.corr.fastq $newname");
			push(@corrected_seqs, "$newname");
		}elsif (-e("./frag_read_corr.corr.fastq")){
			printl("+ WARNING : fragment-only output of ALLPATHS read correction is empty. +\n");
		}else{
			printl("+ WARNING : fragment-only output of ALLPATHS read correction does not exist. +\n");
		}
		
		#PAIRED OUTPUT
		if ( (-s("$expected_pair_1")) && (-s("$expected_pair_2")) ){
			#file found and contains sequence5
			if ($FRAGMENTONLY == 1){
				my $newname = "$DIR/read_correct/" . "$filename" . ".read_corr.frag1.fastq";
				runsys("mv $expected_pair_1 $newname");
				push(@corrected_seqs, "$newname");
				$newname = "$DIR/read_correct/" . "$filename" . ".read_corr.frag2.fastq";
				runsys("mv $expected_pair_2 $newname");
				push(@corrected_seqs, "$newname");
			}
			else{
				my $newout = "$DIR/read_correct/" . "$filename" . "read_corr.pairs.fastq";
				
				#fix qual id line and shuffle
				my $newer_fastq = "$DIR/read_correct/" . "$filename" . ".read_corr.frag1.id.fastq";
				runsys("$AUTOMATON_INSTALL/fastq_ID_fixer.pl -i $expected_pair_1 > $newer_fastq");
				runsys("sed \'/^\$/d\' $newer_fastq > tmp");
				runsys("mv tmp $newer_fastq");
				if (!(-s("$newer_fastq"))){
					printl("Uniqueness check post-processing failed during fastq ID fixing on $newer_fastq .  See logs.\nExiting...\n");
					exit;
				}
				$expected_pair_1 = $newer_fastq;
				$newer_fastq = "$DIR/read_correct/" . "$filename" . ".read_corr.frag2.id.fastq";
				runsys("$AUTOMATON_INSTALL/fastq_ID_fixer.pl -i $expected_pair_2 > $newer_fastq");
				runsys("sed \'/^\$/d\' $newer_fastq > tmp");
				runsys("mv tmp $newer_fastq");
				if (!(-s("$newer_fastq"))){
					printl("Uniqueness check post-processing failed during fastq ID fixing on $newer_fastq .  See logs.\nExiting...\n");
					exit;
				}
				$expected_pair_2 = $newer_fastq;
				runsys("$AUTOMATON_INSTALL/shuffleSequences_fastq.pl $expected_pair_1 $expected_pair_2 $newout");
				if (-s("$newout")){ push(@corrected_seqs, "$newout"); }
				else{ printl("+ ERROR : Failed to interleave pairs.  See logs. +\nExiting...\n"); exit; }
			}
		}elsif ( !(-s("$expected_pair_1")) && (-s("expected_pair_2")) ){
			printl("+ WARNING : No reverse paired-ends exported from ALLPATHS read corrector. +\n");
			my $newname = "$DIR/read_correct/" . "$filename" . ".read_corr.frag1.fastq";
			runsys("mv $expected_pair_1 $newname");
			push(@corrected_seqs, "$newname");
		}elsif ( (-s("$expected_pair_1")) && !(-s("expected_pair_2")) ){
			printl("+ WARNING : No forward paired-ends exported from ALLPATHS read corrector. +\n");
			my $newname = "$DIR/read_correct/" . "$filename" . ".read_corr.frag2.fastq";
			runsys("mv $expected_pair_2 $newname");
			push(@corrected_seqs, "$newname");
		}else{
			printl("+ WARNING : No paired-end output found. +");
		}	
		
		
		#DEBUG
		printl("DEBUG : Cleaning read_correct directory...\n");
		runsys("rm -rf ./data");
		runsys("rm ./in_*");
		runsys("rm ./PE-*");
		runsys("rm ./frag_read_corr.corr.fastq");
		#DEBUG

	}
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return @corrected_seqs;
}



################################################################################################################################################################################################################################################
# uniq_check_frags
################################################################################################################################################################################################################################################
sub uniq_check_frags {
	my $offset = shift;
	my @query_seqs = @_;
	my @uniq_seqs;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/uniq_check"); 
	chdir("$DIR/uniq_check") or die "Could not change directory to $DIR/uniq_check\n"; 
	
	foreach my $seq (@query_seqs){
		my($filename, $directory, $suffix) = fileparse("$seq");
		my $tmp_uniq = "$DIR" . "/uniq_check/" . "$filename" . ".unique.fastq";
		runsys("$FASTQ_UTIL_INSTALL/fastq filter --adjust $offset --unique $seq > $tmp_uniq");
		if (-s("$tmp_uniq")){
			push(@uniq_seqs, "$tmp_uniq");
		}else{
			printl("Uniqueness check failed on sample $filename - see logs for details.\nExiting...\n");
			exit;
		}
	}
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return @uniq_seqs;
}


################################################################################################################################################################################################################################################
# uniq_check_pairs
################################################################################################################################################################################################################################################
sub uniq_check_pairs {
	my $offset = shift;
	my $seq = shift;
	my @uniq_seqs;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/uniq_check"); 
	chdir("$DIR/uniq_check") or die "Could not change directory to $DIR/uniq_check\n"; 
	
	my($filename, $directory, $suffix) = fileparse("$seq");
	my $tmp_uniq = "$DIR" . "/uniq_check/" . "$filename" . ".unique.fastq";
	runsys("$FASTQ_UTIL_INSTALL/fastq filter --adjust $offset --unique $seq > $tmp_uniq");
	if (-s("$tmp_uniq")){
		#push(@uniq_seqs, "$tmp_uniq");
		
		#EVALUATE PAIRS VS FRAGS
		runsys("cat $tmp_uniq | grep \'@\' | grep 'corr' | sort | tr \'@\' \' \' | awk \'{print \$1}\' > ./paired_id_list.txt");
		runsys("$AUTOMATON_INSTALL/extract_fastq_as_mates_frags.pl paired_id_list.txt checkPAIR");
		printl("Now evaluating valid pairs following uniqueness run on paired end sample...\n");
		if ( (-s("./checkPAIR.fragment_ids.txt")) && (-s("./checkPAIR.mate_ids.txt")) ){
			#both outputs found
			#frags first
			my $new_fastq = "$DIR/uniq_check/" . "$filename" . ".unique.fragments.fastq";
			runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $tmp_uniq -name ./checkPAIR.fragment_ids.txt -outfile $new_fastq");
			if (-s($new_fastq)){
				$uniq_seqs[0] = $new_fastq;
			}else{
				printl("+ ERROR during extraction of fragments from unique-check output. +\nExiting...\n");
				exit;
			}
			$new_fastq = "$DIR/uniq_check/" . "$filename" . ".unique.pairs.fastq";
			runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $tmp_uniq -name ./checkPAIR.mate_ids.txt -outfile $new_fastq");
			if (-s($new_fastq)){
				$uniq_seqs[1] = $new_fastq;
			}else{
				printl("+ ERROR during extraction of pairss from unique-check output. +\nExiting...\n");
				exit;
			}
		}elsif (-s("./checkPAIR.fragment_ids.txt")){
			#fragments only
			printl("+ WARNING : No mates exported from uniqueness check on paired-end sample [ $filename ]. +\n");
			my $new_fastq = "$DIR/uniq_check/" . "$filename" . ".unique.fragments.fastq";
			runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $tmp_uniq -name ./checkPAIR.fragment_ids.txt -outfile $new_fastq");
			if (-s($new_fastq)){
				$uniq_seqs[0] = $new_fastq;
			}else{
				printl("+ ERROR during extraction of fragments from unique-check output. +\nExiting...\n");
				exit;
			}
		}elsif (-s("./checkPAIR.mate_ids.txt")){
			#pairs only
			printl("+ WARNING : No fragment-only reads exported from uniqueness check on paired-end sample [ $filename ]. This indicates very low duplicity in your query set. +\n");
			my $new_fastq = "$DIR/uniq_check/" . "$filename" . ".unique.pairs.fastq";
			runsys("$MIRA_INSTALL/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile $tmp_uniq -name ./checkPAIR.mate_ids.txt -outfile $new_fastq");
			if (-s($new_fastq)){
				$uniq_seqs[1] = $new_fastq;
			}else{
				printl("+ ERROR during extraction of pairs from unique-check output. +\nExiting...\n");
				exit;
			}
		}
		else{
			printl("+ ERROR : No output found during mate fixing of uniqueness check output for file [ $filename ]. +\nExiting...\n");
			exit;
		}
		
	}else{
		printl("Uniqueness check failed on sample $filename - see logs for details.\nExiting...\n");
		exit;
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return @uniq_seqs;
}


################################################################################################################################################################################################################################################
# dust_check
################################################################################################################################################################################################################################################
sub dust_check {
	my $type = shift;
	my $offset = shift;
	my $SFFinsert = shift;
	my @query_seqs = @_;
	my @return_locs;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/dust_check"); 
	chdir("$DIR/dust_check") or die "Could not change directory to $DIR/dust_check\n"; 
	
	foreach my $seq (@query_seqs){
		my($filename, $directory, $suffix) = fileparse("$seq");
		my $ztmp = 1;
		my $out_tmp_dust; my $out_tmp_dust2;
		my $new_fastq; my $newer_fastq;
		if ($type eq "sff"){ #SFF
			#CONVERT
			my $libname = "$filename" . "_" . "$ztmp";
			my $storename = "$libname" ."Store";
			my $outfile = "$filename" . ".frg";
			if ( $SFFinsert > 0 ){
				my $SD = sprintf( "%.0f", ($SFFinsert/10) );
				runsys("$CELERA_WGS_INSTALL/bin/sffToCA -libraryname $libname -linker titanium -insertsize $SFFinsert $SD -trim chop -output $outfile $seq");
			}else{
				runsys("$CELERA_WGS_INSTALL/bin/sffToCA -libraryname $libname -trim chop -output $outfile $seq");
			}
			runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $storename $outfile");
			runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfastq $libname $storename");
			#TEST _ DEBUG
			
			
			my $out_tmp =  "$libname" . ".paired.fastq";
			$out_tmp_dust = "$DIR/dust_check/" . "$libname" . ".paired.dust_masked.fastq";
			my $out_tmp2 =  "$libname" . ".unmated.fastq";
			$out_tmp_dust2 = "$DIR/dust_check/" . "$libname" . ".unmated.dust_masked.fastq";

			if ( ( -s("$out_tmp") ) && ( -s("$out_tmp2") ) ){
				#454 PAIRS + FRAGMENTS
				#VERSION 2.0
                           runsys("ln -s $out_tmp $out_tmp_dust");
                           runsys("ln -s $out_tmp2 $out_tmp_dust2");
				
			}
			elsif (-s("$out_tmp")){
				#454 PAIRS ONLY
				#VERSION 2.0
                           runsys("ln -s $out_tmp $out_tmp_dust");
				$out_tmp_dust2 = "";
			}
			elsif (-s("$out_tmp2")){
				#454 FRAGMENTS ONLY
				$out_tmp_dust = "";
				my $new = "$out_tmp2" . ".qv.fastq";
				runsys("quality_trim -n -f 33 -c 0 -r $out_tmp2 -o $new");
				$out_tmp2 = $new;
				my $out_tmp2_idfix = "$DIR/dust_check/" . "$libname" . ".unmated.id.fastq";
				runsys("$AUTOMATON_INSTALL/fastq_ID_fixer.pl -i $out_tmp2 > $out_tmp2_idfix");
				runsys("sed \'/^\$/d\' $out_tmp2_idfix > tmp");
				runsys("mv tmp $out_tmp2_idfix");
				if (!(-s("$out_tmp2_idfix"))){
						printl("DUST post-processing failed during fastq ID fixing.  See logs.\nExiting...\n");
						exit;
				}
				
				my $new_seq = "$out_tmp2" . ".fasta";
				#my $new_qual = "$out_tmp2" . ".dust_masked.fasta.qual";
				printl("EXTRACTING SEQUENCES AND QUALITY BEFORE DUST CHECK\n"); 
				
				#FRAGMENT-ONLY
				if ($out_tmp2_idfix){
					open(GC,"< $out_tmp2_idfix");
					my $count=0;
					open (SEQOUT, '> out.seq');
					open (QUALOUT, '> out.qual');
				    while (defined(my $line = <GC>)){ #get each line
						chomp $line;
						if ($count == 0){
							my $tmp = reverse($line);
							chop($tmp);
							my $tmp2 = reverse($tmp);
							my $line = ">" . "$tmp2";
							print SEQOUT "$line\n";
						}elsif ($count == 1){
							print SEQOUT "$line\n";
						}elsif ($count == 2){ 
							my $tmp = reverse($line);
							chop($tmp);
							my $tmp2 = reverse($tmp);
							my $line = ">" . "$tmp2";
							print QUALOUT "$line\n";
					    }else{
					    	#system("echo \"$line\" >> $new_qual");
					    	print QUALOUT "$line\n";
					    	$count = -1;
					    }
						$count++;
				    }
				    close SEQOUT;
					close QUALOUT;
				}
				runsys("mv out.seq $new_seq");
				my $newer_seq = "$out_tmp2" . ".dust_masked.fasta";
				my $new_qual = "$out_tmp2" . ".dust_masked.qual";
				runsys("mv out.qual $new_qual");
				if ($new_seq && $new_qual){
					runsys("$DUST_INSTALL/dust $new_seq > $newer_seq");
				}else{
					printl("+ ERROR : Extraction of seq vs qual failed +\n");
				}
				if (-s("$newer_seq")){
					runsys("$AUTOMATON_INSTALL/fastqQual2fastq_0.pl $newer_seq");
				}else{
					printl("DUST low complexity run failed during conversion following run.  See logs.\nExiting...\n");
					exit;
				}
				$new_fastq = "$DIR/dust_check/" . "$out_tmp2" . ".dust_masked.fastq";
				if ((-s("$new_fastq"))){
					printl("Expected file found (debug) : $new_fastq\n");
					$out_tmp_dust2 = $new_fastq;
				}else{ printl("454 fragment processing failed following DUST low complexity mask -- see logs.\nExiting...\n"); exit; }
			}
			#DEBUG
			#exit;
			
		}
		else{ #FASTQ
			my $new_seq = "$filename" . ".fasta";
			my $new_qual = "$new_seq" . ".qual"; 

			
			runsys("$RNNOTATOR_INSTALL/bin/fastq_to_fasta_qual $seq $new_seq $new_qual");
			
			$out_tmp_dust = "$DIR/dust_check/" . "$filename" . "dust_masked.fasta";
			my $out_new_qual = "$DIR/dust_check/" . "$filename" . "dust_masked.qual";
			runsys("mv $new_qual $out_new_qual");
			if (-s("$new_seq")){
				runsys("$DUST_INSTALL/dust $new_seq > $out_tmp_dust");
				if (-s("$out_tmp_dust")){

						runsys("$PERL_INSTALL/perl $AUTOMATON_INSTALL/fastaQual2fastq_64_JM.pl $out_tmp_dust");

				}else{
					printl("DUST low complexity run failed during run.  See logs.\nExiting...\n");
					exit;
				}
				$new_fastq = "$DIR/dust_check/" . "$filename" . "dust_masked.fastq";
				if (!(-s("$new_fastq"))){

						printl("DUST post-processing failed during fasta->fastq conversion.  See logs.\nExiting...\n");
						exit;

				}
				#fix id's
				$newer_fastq = "$DIR/dust_check/" . "$filename" . "dust_masked_id.fastq";
				runsys("$AUTOMATON_INSTALL/fastq_ID_fixer.pl -i $new_fastq > $newer_fastq");
				runsys("sed \'/^\$/d\' $newer_fastq > tmp");
				runsys("mv tmp $newer_fastq");
				if (!(-s("$newer_fastq"))){
					printl("DUST post-processing failed during fastq ID fixing.  See logs.\nExiting...\n");
					exit;
				}
				$out_tmp_dust = $newer_fastq;
			}else{
				printl("DUST low complexity run failed due to failed fastq->fasta conversion.  See logs.\nExiting...\n");
				exit;
			}
		}
		$ztmp++;
		
		my $tmp_count2; my $tmp_count3;
		if ($out_tmp_dust){
			push (@return_locs, $out_tmp_dust);
			$tmp_count2 = `$CLC_NGS_INSTALL/sequence_info -r $out_tmp_dust | grep \"Number of N\" | awk \'{print \$4}\'`; 
			chomp $tmp_count2;
		}else{
			$tmp_count2 = 0;
		}
		if ($out_tmp_dust2){
			push (@return_locs, $out_tmp_dust2);
			$tmp_count3 = `$CLC_NGS_INSTALL/sequence_info -r $out_tmp_dust2 | grep \"Number of N\" | awk \'{print \$4}\'`; 
			chomp $tmp_count3;
		}else{
			$tmp_count3 = 0;
		}
		my $tmp_count1 = `$CLC_NGS_INSTALL/sequence_info -r $seq | grep \"Number of N\" | awk \'{print \$4}\'`; chomp $tmp_count1;
		if ($tmp_count3){ $tmp_count2 += $tmp_count3; }
		printl("LOW COMPLEXITY CHECK ON FILE $filename :\n");
		printl("N basecounts before DUST    : $tmp_count1\n");
		printl("N basecounts following DUST : $tmp_count2\n");
		my $tmp_count = $tmp_count2 - $tmp_count1;
		printl("DUST-masked base count      : $tmp_count\n");
		$totalDUSTchanges += $tmp_count;
	}
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return @return_locs;
}

################################################################################################################################################################################################################################################
# qv_trim
################################################################################################################################################################################################################################################
sub qv_trim_frags {
	my $offset = shift;
	my @query_seqs = @_;
	my @trimmed_seqs;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/qv_trim"); 
	chdir("$DIR/qv_trim") or die "Could not change directory to $DIR/qv_trim\n"; 
	
	foreach my $seq (@query_seqs){
		my($filename, $directory, $suffix) = fileparse("$seq");
		my $new_out = "$DIR/qv_trim" . "/" . "$filename" . ".qvTRIM.fragments.fastq";
		my $captured_offset = `$CLC_NGS_INSTALL/quality_trim -n -f $offset -c 18 -r $seq -o $new_out | grep Quality | awk \'{print \$5}\'`;
		chomp $captured_offset;
		runsys("$CLC_NGS_INSTALL/quality_trim -f $offset -n -c 18 -r $seq -o $new_out");
		if (-s("$new_out")){
			push(@trimmed_seqs, "$new_out");
		}else{
			printl("Modifying QV's to match expected range...\n");
			if ($captured_offset < 38){
				my $new_offset = 41 - $captured_offset;
				runsys("$CLC_NGS_INSTALL/quality_trim -f $new_offset -c 18 -r $seq -o $new_out");
				if (-s("$new_out")){
					push(@trimmed_seqs, "$new_out");
				}else{
					printl("Error trimming file $filename -- See Logs...\nExiting...\n\n");
					exit;
				}
			}elsif( ($captured_offset > 42) ){
				my $new_offset = 33;
				runsys("$CLC_NGS_INSTALL/quality_trim -f $new_offset -c 18 -r $seq -o $new_out");
				if (-s("$new_out")){
					push(@trimmed_seqs, "$new_out");
				}else{
					printl("Error trimming file $filename -- See Logs...\nExiting...\n\n");
					exit;
				}
			}else{
				printl("+ WARNING : NO FRAGMENTS KEPT IN OUTPUT +\n");
			}
		}
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	
	return @trimmed_seqs;
}

################################################################################################################################################################################################################################################
# qv_trim
################################################################################################################################################################################################################################################
sub qv_trim_pairs {
	my $offset = shift;
	my $query_seq = shift;
	my @trimmed_seqs;
	
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	runsys("mkdir $DIR/qv_trim"); 
	chdir("$DIR/qv_trim") or die "Could not change directory to $DIR/qv_trim\n"; 
	
	my($filename, $directory, $suffix) = fileparse("$query_seq");
	my $new_frags = "$DIR/qv_trim" . "/" . "$filename" . ".qvTRIM.fragments.fastq";
	my $new_pairs = "$DIR/qv_trim" . "/" . "$filename" . ".qvTRIM.pairs.fastq";
	runsys("$CLC_NGS_INSTALL/quality_trim -n -f $offset -c 18 -r $query_seq -o $new_frags -p $new_pairs");
	if ( (-s("$new_frags")) && (-s("$new_frags")) ){
		$trimmed_seqs[0] = "$new_frags";
		$trimmed_seqs[1] = "$new_pairs";
	}
	elsif (-s("$new_frags")) {
		$trimmed_seqs[0] = "$new_frags";
	}
	elsif (-s("$new_pairs")) {
		$trimmed_seqs[1] = "$new_pairs";
	}
	else{
		printl("Modifying QV's to match expected range...\n");
		my $captured_offset = `$CLC_NGS_INSTALL/quality_trim -f $offset -c 18 -r $query_seq -o $new_frags -p $new_pairs | grep Quality | awk \'{print \$5}\'`;
		chomp $captured_offset;
		
		my $new_offset;
		if ($captured_offset < 38){
			$new_offset = 41 - $captured_offset;
		}elsif( ($captured_offset > 42) ){
			$new_offset = 33;
		}else{
			printl("Acceptable QV offset could not be found for file $filename\nExiting...\n");
			exit;
		}
				
		runsys("$CLC_NGS_INSTALL/quality_trim -n -f $new_offset -c 18 -r $query_seq -o $new_frags -p $new_pairs");
		if ( (-s("$new_frags")) && (-s("$new_frags")) ){
			$trimmed_seqs[0] = "$new_frags";
			$trimmed_seqs[1] = "$new_pairs";
		}
		elsif (-s("$new_frags")) {
			$trimmed_seqs[0] = "$new_frags";
		}
		elsif (-s("$new_pairs")) {
			$trimmed_seqs[1] = "$new_pairs";
		}
		else{
			printl("Error trimming file $filename -- See Logs...\nExiting...\n\n");
			exit;
		}
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	
	return @trimmed_seqs;
}

################################################################################################################################################################################################################################################
# assemble_clc
################################################################################################################################################################################################################################################
sub assemble_clc {
	my $SFFinsert = shift; #1=presub #0=contam
	my $ILLUMinsert = shift; #sff fasta fastq
	my $prefix = shift;
	my %dataset_hash = @_;
	
	#EXTRACT HASH
	my @tmp1; my @tmp2; my @tmp3; my @tmp4; my @tmp5;
	my @c1; my @c2; my @c3; my @c4; my @c5;
	if ( @{$dataset_hash{0}}[0] ){ @c1 = @{$dataset_hash{0}};	}
	if ( @{$dataset_hash{1}}[0] ){ @c2 = @{$dataset_hash{1}};	}
	if ( @{$dataset_hash{2}}[0] ){@c3 = @{$dataset_hash{2}}; 	}
	if ( @{$dataset_hash{3}}[0] ){ @c4 = @{$dataset_hash{3}};	}
	if ( @{$dataset_hash{4}}[0] ){ @c5 = @{$dataset_hash{4}};	}
	
	
	#SET UP DIRECTORY
	if (!(-e("$DIR/CLC_ASSEMBLIES"))){ mkdir ("$DIR/CLC_ASSEMBLIES"); }
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	chdir("$DIR/CLC_ASSEMBLIES") or die "Could not change directory to $DIR/CLC_ASSEMBLIES\n";
	if (!(-e "./license.properties")){ runsys("ln -s $CLC_NGS_INSTALL/license.properties ."); }
	
	#PROCESS 454 PAIRS
	my @new_sff_pairs;
	my @new_sff_fragments;
	#foreach tmp2
	foreach(@c2){
		my $tmp_p = "$_" . ".pairs.fasta";
		my $tmp_f = "$_" . ".fragments.fasta";
		runsys("$CLC_NGS_INSTALL/split_sequences -i $_ -s $tmp_f -p $tmp_p -d ti");
		if (-s("$tmp_f")){ push(@new_sff_fragments, "$tmp_f"); }
		if (-s("$tmp_p")){ push(@new_sff_pairs, "$tmp_p"); }
	}
	
	#BUILD RUN COMMAND
	my ($sffSD, $sff_min, $sff_max, $illumSD, $illum_min, $illum_max);
	my $tmp_out = "$prefix" . "_clc_contigs.fasta";
	if ($SFFinsert){
		$sffSD = sprintf( "%.0f", ($SFFinsert/10)) ;
		$sff_min = $SFFinsert - $sffSD;
		$sff_max = $SFFinsert + $sffSD;
	}
	if ($ILLUMinsert){
		$illumSD = sprintf( "%.0f", ($ILLUMinsert/10)) ;
		$illum_min = $ILLUMinsert - $illumSD;
		$illum_max = $ILLUMinsert + $illumSD;
	}
	
	#COMMAND
	my $assemb_cmd;
	if ( (@c1 > 0) && (@new_sff_pairs > 0) && (@new_sff_fragments > 0) && (@c3 > 0) && (@c4 > 0) ){
		$assemb_cmd = "$CLC_NGS_INSTALL/clc_novo_assemble --cpus 4 -o $tmp_out -q @c1 @new_sff_fragments @c3 -p fb ss $sff_min $sff_max @new_sff_pairs -p fb ss $illum_min $illum_max @c4";
	}
	else{
		$assemb_cmd = "$CLC_NGS_INSTALL/clc_novo_assemble --cpus 4 -o $tmp_out -q ";
		if (@c1 > 0){ $assemb_cmd .= "@c1 "; }
		if (@new_sff_fragments > 0){ $assemb_cmd .= "@new_sff_fragments "; }
		if (@c3 > 0){ $assemb_cmd .= "@c3 "; }
		if (@c5 > 0){ $assemb_cmd .= "@c5 "; }
		if (@new_sff_pairs > 0){ $assemb_cmd .= "-p fb ss $sff_min $sff_max @new_sff_pairs " }
		if (@c4 > 0){ $assemb_cmd .= "-p fb ss $illum_min $illum_max @c4" }
	}
	
	#EXECUTE
	runsys("$assemb_cmd");
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	
	#LEAVE SUBROUTINE
	my $clc_output = "$DIR/CLC_ASSEMBLIES" . "/" . "$tmp_out";
	printl("Assembled CLC output : $clc_output");
	printl("\n\n");
	return $clc_output; 
}


################################################################################################################################################################################################################################################
# assemble_velvet
################################################################################################################################################################################################################################################
sub assemble_velvet {
	my $SFFinsert = shift;
	my $ILLUMinsert = shift;
	my $format = shift; #sff fasta fastq
	my $prefix = shift;
	my $exp_cov = shift;
	my %dataset_hash = @_;

	#SET UP DIRECTORY
	if (!(-e("$DIR/VELVET_ASSEMBLIES"))){ mkdir ("$DIR/VELVET_ASSEMBLIES"); }
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	chdir("$DIR/VELVET_ASSEMBLIES") or die "Could not change directory to $DIR/VELVET_ASSEMBLIES\n";
	
	#PROCESS 454 PAIRS
	my @new_sff_pairs;
	my @new_sff_fragments;
	if ($format == "sff"){
		if ( @{$dataset_hash{1}} ){
			foreach( @{$dataset_hash{1}} ){
				my $tmp_p = "$_" . ".pairs.fasta";
				my $tmp_f = "$_" . ".fragments.fasta";
				runsys("$CLC_NGS_INSTALL/split_sequences -i $_ -s $tmp_f -p $tmp_p -d ti");
				if (-s("$tmp_f")){ push(@new_sff_fragments, "$tmp_f"); }
				if (-s("$tmp_p")){ push(@new_sff_pairs, "$tmp_p"); }
				if ( !(-s("$tmp_f")) && !(-s("$tmp_p")) ){
					printl("No split 454 output found during paired end prep for Velvet assembly\n( file : $_ )\nExiting...");
					exit;
				}
			}
		}
	}

#ASSEMB SWITCH	
if ($no_assemb_switch == 0){
	
	#VELVETH INITIAL STEP!!!!
	my $paired_count = 0;
	my @paired_insert = ();
	my $velvetdir = "$DIR/VELVET_ASSEMBLIES/" . "$prefix" . "/";
	my $velvetH_cmd = "$VELVET_INSTALL/velveth $velvetdir 25";
	if ( @{$dataset_hash{0}}[0] ){ 
		#454 FRAGMENTS
		$velvetH_cmd .= " -fastq -long @{$dataset_hash{0}}";
	}
	if ( @{$dataset_hash{2}}[0] ){
		#ILLUMINA FRAGMENTS
		if ( @{$dataset_hash{0}}[0] ){ 
			#IF 454 FRAGS ALSO EXIST (no flags required)
				$velvetH_cmd .= " @{$dataset_hash{2}}";
		}else{
			#IF NOT
				$velvetH_cmd .= " -fastq -long @{$dataset_hash{2}}";
		}
	}	
	if ( @{$dataset_hash{4}}[0] ){ 
		#fasta files
		$velvetH_cmd .= " -fasta -long @{$dataset_hash{4}}";
	}
	if ( @{$dataset_hash{1}}[0] ){ 
		#454 PAIRS
		if ($format == "sff"){
			if ($new_sff_fragments[0]){
				$velvetH_cmd .= " -fasta -long @new_sff_fragments";
			}
			if ($new_sff_pairs[0]){
				my $pair_string;
				if ( $paired_count == 0 ){ $pair_string = "-longPaired"; }
				else{ $pair_string = "-longPaired" . "$paired_count"; }
				$velvetH_cmd .= " -fasta $pair_string @new_sff_pairs";
				$paired_insert[$paired_count] = $SFFinsert;
				$paired_count++;
			}
		}else{
			$velvetH_cmd .= " -fastq -long @{$dataset_hash{1}}";
		}
	}
	if ( @{$dataset_hash{3}}[0] ){ 
		#ILLUMINA PAIRS
		my $pair_string;
		if ( $paired_count == 0 ){ $pair_string = "-longPaired"; }
		else{ $pair_string = "-longPaired" . "$paired_count"; }
		$velvetH_cmd .= " -fastq $pair_string @{$dataset_hash{3}}";
		$paired_insert[$paired_count] = $ILLUMinsert;
		$paired_count++;
	}
	#RUN VELVETH
	printl("\n\nRUNNING VELVET ASSEMBLY PT. 1 (velveth) :\n");
	runsys("$velvetH_cmd");
	
	

	#VELVETG FINAL STEP
	my $velvetG_cmd = "$VELVET_INSTALL/velvetg $prefix -exp_cov $exp_cov";
	my $r = 0;
	if (-e("$velvetdir")){
		foreach(@paired_insert){
			my $sd;
			if ($r == 0){
				$sd = sprintf( "%.0f", ($_/10) );
				$velvetG_cmd .= " -ins_length_long $_ -ins_length_long_sd $sd";
			}else{
				$sd = sprintf( "%.0f", ($_/10) );
				my $str1 = "-ins_length_long" . "$r";
				my $str2 = "-ins_length_long_sd" . "$r";
				$velvetG_cmd .= " $str1 $_ $str2 $sd";
			}
			$velvetG_cmd .= " -amos_file yes -unused_reads yes";
			$r++;
		}
		#RUN VELVETG
		printl("\n\nRUNNING VELVET ASSEMBLY PT. 2 (velvetg) :\n");
		runsys("$velvetG_cmd");
	}else{
		printl("+ ERROR : VELVETH Run Failed! No output directory found at: +\n $velvetdir\n\nExiting...");
		exit;
	}

	#CHECK FOR OUTPUT
	my $assembly_loc = "$DIR/VELVET_ASSEMBLIES/" . "$prefix" . "/contigs.fa";
	if (!(-s("$assembly_loc"))){
		printl("+ ERROR : No velvet output found! +\nExiting...\n");
		#exit;
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return $assembly_loc;
#END ASSEMB SWITCH	
}else{ return "NONE"; }
}



################################################################################################################################################################################################################################################
# assemble_newbler
################################################################################################################################################################################################################################################
sub assemble_newbler {
	my $SFFinsert = shift;
	my $ILLUMinsert = shift;
	my $format = shift; #sff fasta fastq
	my $prefix = shift;
	my $offset = shift;
	my %dataset_hash = @_;

	#SET UP DIRECTORY
	if (!(-e("$DIR/NEWBLER_ASSEMBLIES"))){ mkdir ("$DIR/NEWBLER_ASSEMBLIES"); }
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	chdir("$DIR/NEWBLER_ASSEMBLIES") or die "Could not change directory to $DIR/NEWBLER_ASSEMBLIES\n";
	
	printl("\n\nPREPARING FOR NEWBLER ASSEMBLY:\n");
	
	#BUILD RUN COMMAND
	my $newbler_run_cmd = "$SEQ454_INSTALL/runAssembly -o $DIR/NEWBLER_ASSEMBLIES/assembly -force -p";
	my @new_query_seqs;	
		
	my $mock_count = 0;
	
		#454 FILES FIRST
		if ($format eq "sff"){
			#SFF INPUT FILES
			if ( @{$dataset_hash{0}} ){
				#FRAGMENTS
				$newbler_run_cmd .= " @{$dataset_hash{0}}";
			}
			if ( @{$dataset_hash{1}} ){
				#PAIRS
				$newbler_run_cmd .= " @{$dataset_hash{1}}";
			}
		}else{
			#FASTQ INPUT FILES
			if ( @{$dataset_hash{0}} ){
				#FRAGMENTS
				foreach ( @{$dataset_hash{0}} ){
					printl("Converting 454 (fastq) to 454 (fna+qual) for Newbler assembly:\n");
					#CONVERT TO FNA+QUAL
					my $libname = "frag" . "$mock_count";
					my $mock_out = "$libname" . ".mock.frg";
					my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
					if ($offset == 64){
						$convert_cmd .= " -type illumina";
					}elsif($offset == 33){
						$convert_cmd .= " -type sanger";
					}elsif($offset == 0){
						$convert_cmd .= " -type solexa";
					}else{
						$convert_cmd .= " -type illumina";	
					}
					$convert_cmd .= " -technology 454";
					$convert_cmd .= " -reads $_ > $mock_out";
					runsys("$convert_cmd");
					if (-s("$mock_out")){
						my $libstore = "$libname" . "Store";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
						if (-e("$libstore")){
							runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpnewbler $libname $libstore");
							my $new_out = "$libname" . ".fna";
							if (-s($new_out)){
								push(@new_query_seqs, "$new_out");
							}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
					$mock_count++;
				}
			}
			if ( @{$dataset_hash{1}} ){
				#PAIRS
				#(VERSION 2.0)
				foreach ( @{$dataset_hash{1}} ){
					printl("Converting 454 (fastq) to 454 (fna+qual) for Newbler assembly:\n");
					#CONVERT TO FNA+QUAL
					my $libname = "pair" . "$mock_count";
					my $mock_out = "$libname" . ".mock.frg";
					my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
					#$convert_cmd .= " -type illumina";
					if ($offset == 64){
						$convert_cmd .= " -type illumina";
					}elsif($offset == 33){
						$convert_cmd .= " -type sanger";
					}elsif($offset == 0){
						$convert_cmd .= " -type solexa";
					}else{
						$convert_cmd .= " -type illumina";	
					}
					$convert_cmd .= " -technology 454";
					my $SFFsd = sprintf( "%.0f", ($SFFinsert/10));
					$convert_cmd .= " -insertsize $SFFinsert $SFFsd";
					$convert_cmd .= " -mates $_ > $mock_out";

					runsys("$convert_cmd");
					if (-s("$mock_out")){
						my $libstore = "$libname" . "Store";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
						if (-e("$mock_out")){
							runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpnewbler $libname $libstore");
							my $new_out = "$libname" . ".fna";
							if (-s($new_out)){
								push(@new_query_seqs, "$new_out");
							}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
					$mock_count++;
				}
			}
		}
		
		#ILLUMINA FILES SECOND
		
		if ( @{$dataset_hash{2}} ){	
			foreach ( @{$dataset_hash{2}} ){
			#ILLUMINA FRAGMENTS
				printl("Converting ILLUMINA (fastq) to 454 (fna+qual) for Newbler assembly:\n");
				#CONVERT TO FNA+QUAL
				my $libname = "frag" . "$mock_count";
				my $mock_out = "$libname" . ".mock.frg";
				my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
				$convert_cmd .= " -type illumina";
				$convert_cmd .= " -technology illumina";
				$convert_cmd .= " -reads $_ > $mock_out";
				runsys("$convert_cmd");
				if (-s("$mock_out")){
					my $libstore = "$libname" . "Store";
					runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
					if (-e("$libstore")){
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpnewbler $libname $libstore");
						my $new_out = "$libname" . ".fna";
						if (-s($new_out)){
							push(@new_query_seqs, "$new_out");
						}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
				}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
				$mock_count++;
			}
		}
		if ( @{$dataset_hash{3}} ){
			foreach ( @{$dataset_hash{3}} ){
			#ILLUMINA PAIRS
				printl("Converting 454 (fastq) to 454 (fna+qual) for Newbler assembly:\n");
				#CONVERT TO FNA+QUAL
				my $libname = "pair" . "$mock_count";
				my $mock_out = "$libname" . ".mock.frg";
				my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
				$convert_cmd .= " -type illumina";
				$convert_cmd .= " -technology 454";
				my $ILLUMsd = sprintf( "%.0f", ($ILLUMinsert/10));
				$convert_cmd .= " -insertsize $ILLUMinsert $ILLUMsd";
				$convert_cmd .= " -mates $_ > $mock_out";
				runsys("$convert_cmd");
				if (-s("$mock_out")){
					my $libstore = "$libname" . "Store";
					runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
					if (-e("$libstore")){
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpnewbler $libname $libstore");
						my $new_out = "$libname" . ".fna";
						if (-s($new_out)){
							push(@new_query_seqs, "$new_out");
						}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
				}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
				$mock_count++;
			}
		}
		if ( @{$dataset_hash{4}} ){ #fragment fasta files
			foreach( @{$dataset_hash{4}} ){
				push(@new_query_seqs, "$_");
			}
		}
		
#ASSEMB SWITCH	
if ($no_assemb_switch == 0){
	
	#ALL FILES HANDLED - LAUNCH ASSEMB	
	$newbler_run_cmd .= " @new_query_seqs";
	printl("Launching Newbler assembly...\n");
	runsys("$newbler_run_cmd");
	
	#FIND OUTPUT
	my $assembly_loc = "$DIR/NEWBLER_ASSEMBLIES/assembly/454AllContigs.fna";
	if (!(-s("$assembly_loc"))){
		printl("Newbler output not found in expected location:\n$assembly_loc\nSee logs.  Exiting...\n");
		#exit;
		$assembly_loc = "FAILED";
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return $assembly_loc;
	
#END ASSEMB SWITCH	
}else{ return "NONE"; }

}


################################################################################################################################################################################################################################################
# assemble_celera
################################################################################################################################################################################################################################################
sub assemble_celera {
	my $SFFinsert = shift;
	my $ILLUMinsert = shift;
	my $format = shift; #sff format : sff fasta fastq
	my $prefix = shift;
	my $offset = shift;
	my %dataset_hash = @_;

	#SET UP DIRECTORY
	if (!(-e("$DIR/CELERA_ASSEMBLIES"))){ mkdir ("$DIR/CELERA_ASSEMBLIES"); }
	my $cur_pwd = `pwd`; chomp $cur_pwd;
	chdir("$DIR/CELERA_ASSEMBLIES") or die "Could not change directory to $DIR/CELERA_ASSEMBLIES\n";
	
	printl("\n\nPREPARING FOR CELERA ASSEMBLY:\n");
	
	#PREP INPUTFILES!
	my $mock_count = 0;
	my @new_query_seqs = ();
	
	#454 FILES FIRST
		if ($format eq "sff"){
			#SFF INPUT FILES
			if ( @{$dataset_hash{0}} ){
				#FRAGMENTS
				foreach( @{$dataset_hash{0}} ){
					my $libname = "frag" . "$mock_count";
					my $mock_out = "$libname" . ".mock.frg";
					my $convert_cmd = "$CELERA_WGS_INSTALL/bin/sffToCA";
					$convert_cmd .= " -trim chop -libraryname $libname -output $mock_out $_";
					runsys("$convert_cmd");
					if (-s("$mock_out")){
						my $libstore = "$libname" . "Store";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
						if (-e("$libstore")){
							my $new_out = "$DIR/CELERA_ASSEMBLIES/" . "$libname" . ".frg";
							runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfrg $libstore > $new_out");
							if (-s($new_out)){
								push(@new_query_seqs, "$new_out");
							}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
					$mock_count++;
				}
			}
			if ( @{$dataset_hash{1}} ){
				#PAIRS
				foreach( @{$dataset_hash{1}} ){
					my $libname = "frag" . "$mock_count";
					my $mock_out = "$libname" . ".mock.frg";
					my $convert_cmd = "$CELERA_WGS_INSTALL/bin/sffToCA";
					my $SFFsd = sprintf( "%.0f", ($SFFinsert/10));
					$convert_cmd .= " -linker titanium -insertsize $SFFinsert $SFFsd";
					$convert_cmd .= " -trim chop -libraryname $libname -output $mock_out $_";
					runsys("$convert_cmd");
					if (-s("$mock_out")){
						my $libstore = "$libname" . "Store";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
						if (-e("$libstore")){
							my $new_out = "$DIR/CELERA_ASSEMBLIES/" . "$libname" . ".frg";
							runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfrg $libstore > $new_out");
							if (-s($new_out)){
								push(@new_query_seqs, "$new_out");
							}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
					$mock_count++;
				}
			}
		}else{
			#FASTQ INPUT FILES
			if ( @{$dataset_hash{0}} ){
				foreach ( @{$dataset_hash{0}} ){
					printl("Converting 454 (fastq) to FRG for Celera assembly:\n");
					#CONVERT TO FNA+QUAL
					my $libname = "frag" . "$mock_count";
					my $mock_out = "$libname" . ".mock.frg";
					my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
					#$convert_cmd .= " -type illumina";
					if ($offset == 64){
						$convert_cmd .= " -type illumina";
					}elsif($offset == 33){
						$convert_cmd .= " -type sanger";
					}elsif($offset == 0){
						$convert_cmd .= " -type solexa";
					}else{
						$convert_cmd .= " -type illumina";	
					}
					$convert_cmd .= " -technology 454";
					$convert_cmd .= " -reads $_ > $mock_out";
					runsys("$convert_cmd");
					if (-s("$mock_out")){
						my $libstore = "$libname" . "Store";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
						if (-e("$libstore")){
							my $new_out = "$DIR/CELERA_ASSEMBLIES/" . "$libname" . ".frg";
							runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfrg $libstore > $new_out");
							if (-s($new_out)){
								push(@new_query_seqs, "$new_out");
							}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
					$mock_count++;
				}
			}
			if ( @{$dataset_hash{1}} ){
				foreach ( @{$dataset_hash{1}} ){
					printl("Converting 454 (fastq) to FRG for Celera assembly:\n");
					#CONVERT TO FNA+QUAL
					my $libname = "pair" . "$mock_count";
					my $mock_out = "$libname" . ".mock.frg";
					my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
					#$convert_cmd .= " -type illumina";
					if ($offset == 64){
						$convert_cmd .= " -type illumina";
					}elsif($offset == 33){
						$convert_cmd .= " -type sanger";
					}elsif($offset == 0){
						$convert_cmd .= " -type solexa";
					}else{
						$convert_cmd .= " -type illumina";	
					}
					$convert_cmd .= " -technology 454";
					my $SFFsd = sprintf( "%.0f", ($SFFinsert/10));
					$convert_cmd .= " -insertsize $SFFinsert $SFFsd";
					$convert_cmd .= " -mates $_ > $mock_out";
					runsys("$convert_cmd");
					if (-s("$mock_out")){
						my $libstore = "$libname" . "Store";
						runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
						if (-e("$mock_out")){
							my $new_out = "$DIR/CELERA_ASSEMBLIES/" . "$libname" . ".frg";
							runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfrg $libstore > $new_out");
							if (-s($new_out)){
								push(@new_query_seqs, "$new_out");
							}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
						}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
					}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
					$mock_count++;
				}
			}
		}
		
	#NOW ILLUMINA INPUT FILES
	if ( @{$dataset_hash{2}} ){	
		foreach ( @{$dataset_hash{2}} ){
		#ILLUMINA FRAGMENTS
			printl("Converting ILLUMINA (fastq) to FRG for Celera assembly:\n");
			#CONVERT TO FNA+QUAL
			my $libname = "frag" . "$mock_count";
			my $mock_out = "$libname" . ".mock.frg";
			my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
			$convert_cmd .= " -type illumina";
			$convert_cmd .= " -technology illumina";
			$convert_cmd .= " -reads $_ > $mock_out";
			runsys("$convert_cmd");
			if (-s("$mock_out")){
				my $libstore = "$libname" . "Store";
				runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
				if (-e("$libstore")){
					my $new_out = "$DIR/CELERA_ASSEMBLIES/" . "$libname" . ".frg";
					runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfrg $libstore > $new_out");
					if (-s($new_out)){
						push(@new_query_seqs, "$new_out");
					}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
				}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
			}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
			$mock_count++;
		}
	}
	if ( @{$dataset_hash{3}} ){
		foreach ( @{$dataset_hash{3}} ){
		#ILLUMINA PAIRS
			printl("Converting 454 (fastq) to 454 (fna+qual) for Newbler assembly:\n");
			#CONVERT TO FNA+QUAL
			my $libname = "pair" . "$mock_count";
			my $mock_out = "$libname" . ".mock.frg";
			my $convert_cmd = "$CELERA_WGS_INSTALL/bin/fastqToCA -libraryname $libname";
			$convert_cmd .= " -type illumina";
			$convert_cmd .= " -technology 454";
			my $ILLUMsd = sprintf( "%.0f", ($ILLUMinsert/10));
			$convert_cmd .= " -insertsize $ILLUMinsert $ILLUMsd";
			$convert_cmd .= " -mates $_ > $mock_out";
			runsys("$convert_cmd");
			if (-s("$mock_out")){
				my $libstore = "$libname" . "Store";
				runsys("$CELERA_WGS_INSTALL/bin/gatekeeper --frgminlen 40 -o $libstore $mock_out");
				if (-e("$libstore")){
					my $new_out = "$DIR/CELERA_ASSEMBLIES/" . "$libname" . ".frg";
					runsys("$CELERA_WGS_INSTALL/bin/gatekeeper -dumpfrg $libstore > $new_out");
					if (-s($new_out)){
						push(@new_query_seqs, "$new_out");
					}else{ printl("Gatekeeper dump failed -- see logs.\nExiting...\n"); exit; }
				}else{ printl("Gatekeeper build failed -- see logs.\nExiting...\n"); exit; }
			}else{ printl("Gatekeeper mock fastq build failed -- see logs.\nExiting...\n"); exit; }
			$mock_count++;
		}
	}

#ASSEMB SWITCH	
if ($no_assemb_switch == 0){
	
	#BUILD ASSEMBLY COMMAND AND RUN
	printl("Running assembly:\n");
	runsys("echo \'doToggle=1\' > ./this.spec");
	runsys("echo \'createAGP=1\' >> ./this.spec");
	my $celera_run_cmd = "$CELERA_WGS_INSTALL/bin/runCA -d $DIR/CELERA_ASSEMBLIES/assembly -p $prefix -s $DIR/CELERA_ASSEMBLIES/this.spec @new_query_seqs";
	runsys("$celera_run_cmd");
	
	#FIND OUTPUT
	my $assembly_loc = "$DIR/CELERA_ASSEMBLIES/assembly/10-toggledAsm/9-terminator/" . "$prefix" . ".ctg.fasta";
	if (!(-s("$assembly_loc"))){
		printl("Toggled Celera output not found.  Looking for non-toggled...\n");
		$assembly_loc = "$DIR/CELERA_ASSEMBLIES/assembly/9-terminator/" . "$prefix" . ".ctg.fasta";
		if (!(-s("$assembly_loc"))){
			printl("+ ERROR : No Celera output files found!  See logs.\nExiting...\n");
			#exit;
			$assembly_loc = "FAILED";
		}
	}
	
	chdir("$cur_pwd") or die "Could not change back to directory to $cur_pwd\n";
	return $assembly_loc;
#END ASSEMB SWITCH	
}else{ return "NONE"; }
}



################################################################################################################################################################################################################################################
# MAIN : 
################################################################################################################################################################################################################################################
MAIN : {
	my ($PREFIX, $GENOMESIZE, $BOOLcovcheck, $BOOLuniqcheck, $BOOLdustcheck, $readCOR, $BOOLchaffcheck, $offset, $presub_contam, $contam, @unpaired_454_inputs, @paired_454_inputs, @unpaired_illum_inputs, @paired_illum_inputs, @longmate_illum_inputs, $SFFinsert, $ILLUMinsert);
	$BOOLcovcheck = 0; $BOOLuniqcheck = 0; $BOOLdustcheck = 0; $BOOLchaffcheck = 0;
	$no_assemb_switch = 0;
	my $status = GetOptions(\%Opts, "help!", 'o=s'=> \$DIR, 'g=s'=> \$GENOMESIZE, 'p=s'=> \$PREFIX, 'f=s'=> \$offset, 'precontam=s'=> \$presub_contam, 'contam=s'=> \$contam, "c!", "h!", "d!", "z!", "n!", 'r=s'=> \$readCOR, "u!", 'adapter=s@'=> \@adapters, 'SFFfragment=s@'=> \@unpaired_454_inputs, 'SFFpair=s@'=> \@paired_454_inputs, 'ILLUMfragment=s@'=> \@unpaired_illum_inputs, 'ILLUMpair=s@'=> \@paired_illum_inputs, 'SFFinsert=s'=> \$SFFinsert, 'ILLUMinsert=s'=> \$ILLUMinsert);
	
	#PARSE INPUT------------------------------------------------------------------------------------------------------------------START
	#print %Opts;
	print "\nUSER INPUT:\n";
	if ( exists $Opts{help}){ printl("Help requested:\n"); outhelp(); }	
	if ( exists $Opts{h}){ printl("Help requested:\n"); outhelp(); }
	if ( exists $Opts{d}){ printl("DUST low complexity check requested (with flag -d)\n"); $BOOLdustcheck = 1; }
	if ( exists $Opts{c}){ printl("Coverage Check Requested (with flag -c)\n"); $BOOLcovcheck = 1; }
	if ( exists $Opts{u}){ printl("Uniqueness check requested (with flag -u)\n"); $BOOLuniqcheck = 1; }
	if ( exists $Opts{z}){ printl("Post assembly chaff analysis requested (with flag -z)\n"); $BOOLchaffcheck = 1; }
	if ( exists $Opts{n}){ printl("Assembly following processing disabled (wih flag -n)\n"); $no_assemb_switch = 1; }
	if ( $DIR ){
		if (-e($DIR)){ printl("Output directory : $DIR\n"); }
		else{
			system("mkdir $DIR");
			if (-e($DIR)){ printl("Output directory (newly created) : $DIR\n"); }
			else{
				printl("+ Bad input directory (-o) does not exist or cannot be created. +\n");
				outhelp();
			}
		}
	}else{
		printl("Output directory - no entry, using cwd : $pwd\n");
		$DIR = $pwd;
	}
	if ( $readCOR ){
		if ( $readCOR eq "quake" ){
			printl("Read Correction Requested (with flag -r) Using QUAKE\n");
		}
		elsif( $readCOR eq "allpaths" ){
			printl("Read Correction Requested (with flag -r) Using ALLPATHS\n");
			my $ztmp = `echo \$PYTHONPATH`; chomp $ztmp;
			if (!($ztmp)){
				printl("\n\n+ ERROR : ALLPATHS ENVIRONMENT VARS NOT SET.  PLEASE SET THE FOLLOWING:\n\n");
				printl("CENTOS5\nsetenv PYTHONPATH /usr/local/devel/BCIS/external_software/ALLPATH_EC\n");
				printl("setenv PATH \"/usr/local/devel/BCIS/external_software/ALLPATH_EC/bin:\$PATH\"\n");
				printl("setenv LD_LIBRARY_PATH /usr/local/packages/atlas/lib\n");
				printl("setenv LD_LIBRARY_PATH \"/usr/local/packages/gcc-4.6.2/lib64:\$LD_LIBRARY_PATH\"\n\n");
				printl("CENTOS6\nsetenv PYTHONPATH /usr/local/devel/BCIS/external_software/ALLPATH_EC\n");
				printl("setenv PATH \"/usr/local/devel/BCIS/external_software/ALLPATH_EC/bin:\$PATH\"\n");
				printl("setenv PATH /usr/local/packages/gcc-4.7.1/bin:/usr/local/packages/python-2.7.3/bin:\$PATH\n");
				printl("setenv setenv LD_LIBRARY_PATH /usr/local/packages/gcc-4.7.1/lib64:/usr/local/packages/python-2.7.3/lib\n\n");
				printl("Exiting...\n");
				exit;
			}
		}
		else{
			printl("+ Bad input for read correction method : $readCOR +\nMust be \"quake\" or \"allpaths\"\n");
			outhelp();
		}
	}
	if ( $GENOMESIZE ){
		if ($GENOMESIZE ne "auto"){ 
			if ($GENOMESIZE > 0){
				printl("Predicted Genome size = $GENOMESIZE\n"); 
			}else{
				printl("Incorrect data entry for flag -g. Exiting.\n");
				exit;
			}
		}
		else{ printl("Predicted Genome size will be automatically calculated.\n"); }
	}else{
		printl("+ Bad input for predicted genome size.  Choose \"-g auto\" to use consensus length for N50 calculation. +\n");
		outhelp();
	}
	if ( $PREFIX ){
		printl("Prefix = $PREFIX\n");
	}else{
		printl("Prefix - no entry, using default : OUT\n");
		$PREFIX = "OUT";
	}
	if ( $presub_contam ){
		if (-s("$presub_contam")){
			printl("Presubmission Contaminant Check Reference = $presub_contam\n");
		}else{
			printl(" + Presubmission Contaminant Check Reference does not exist : $presub_contam\n");
			outhelp();
		}
	}else{
		printl("+ No input chosen for optional presubmission contamination check.  This process will not be run. +\n");
	}
	if ( $contam ){
		if (-s("$contam")){ 
			printl("Contaminant Check Reference = $contam\n");
		}else{
			printl(" + Contaminant Check Reference does not exist : $contam\n");
			outhelp();
		}
	}else{
		printl("+ No input chosen for contamination check.  This process will not be run. +\n");
	}
	if ( $offset ){
		printl("Using user input offset = $offset\n");
	}else{
		printl("+ No illumina QV offset chosen with flag (-f).  Using default offset = 0. +\n");
		$offset = 0;
	}
	if ( @adapters ){
		printl("****\nAdapter sequences read:\n");
		foreach (@adapters){
			printl("$_\n");
		}
		printl("--DONE--\n");
	}
	#INPUT summary prep
	my $c1 = 0; my $c2 = 0; my $c3 = 0; my $c4 = 0;
	printl("\n");
	if ( @unpaired_454_inputs ){
		printl("****\nUn-paired 454 Input files read:\n");
		foreach (@unpaired_454_inputs){
			printl("$_\n");
			$c1++;
		}
		print "DONE - $c1 fragment-only 454 input files found\n****\n";
	}else{
		printl("No fragment-only 454 inputs found.\n")
	}
	if ( @paired_454_inputs ){
		printl("****\nPaired 454 Input files read:\n");
		foreach (@paired_454_inputs){
			printl("$_\n");
			$c2++;
		}
		print "DONE - $c2 paired input 454 files found\n****\n";
	}else{
		printl("No paired 454 inputs found.\n")
	}
	if ( @unpaired_illum_inputs ){
		printl("****\nUn-paired Illum Input files read:\n");
		foreach (@unpaired_illum_inputs){
			printl("$_\n");
			$c3++;
		}
		print "DONE - $c3 fragment-only illumina input files found\n****\n";
	}else{
		printl("No fragment-only illumina inputs found.\n")
	}
	if ( @paired_illum_inputs ){
		printl("****\nPaired Illum Input files read:\n");
		foreach (@paired_illum_inputs){
			printl("$_\n");
			$c4++;
		}
		print "DONE - $c4 paired input illumina files found\n****\n";
	}else{
		printl("No paired illumina inputs found.\n")
	}
	if ( $SFFinsert ){
		printl("454 insert size = $SFFinsert\n");
	}elsif ($c2 > 0){
		printl("+ No input given for flag (-SFFinsert) despite supplying 454 pairs. +\n");
		outhelp();
	}else{
		#...No action.
	}
	if ( $ILLUMinsert ){
		printl("Illumina insert size = $ILLUMinsert\n");
	}elsif($c4 > 0){
		printl("+ No input given for flag (-ILLUMinsert) despite supplying illumina pairs. +\n");
		outhelp();
	}else{
		#...No action.
	}
	
	printl("\n");
	#PARSE INPUT-------------------------------------------------------------------------------------------------------------------STOP
	
	#...
	
	#SUMMARIZE INPUT--------------------------------------------------------------------------------------------------------------START
	my %dataset_hash;
	@{$dataset_hash{0}} = @unpaired_454_inputs;
	@{$dataset_hash{1}} = @paired_454_inputs;
	@{$dataset_hash{2}} = @unpaired_illum_inputs;
	@{$dataset_hash{3}} = @paired_illum_inputs;
	@{$dataset_hash{4}} = ();
	printl("\n\n==== INITIAL ALL READ STATS ====\n");
	summarize_reads(%dataset_hash);
	
	my %assembly_hash;
	my %assembly_data;
	my $assembly_hash_count = 0;
	#SUMMARIZE INPUT---------------------------------------------------------------------------------------------------------------STOP	
	
	
	#FILTERING PIPELINE-----------------------------------------------------------------------------------------------------------START
	my $tmp_clc_output;
	if (-e("$DIR")){
		chdir("$DIR") or die "Could not change directory to $DIR\n";
		printl("\n\n==== INITIAL ALL READS ASSEMBLY (CLC) ====\n");
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "001_all_reads", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("001_all_reads", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
		
		if ($GENOMESIZE eq "auto"){
			my $GENOMESIZE = `$CLC_NGS_INSTALL/sequence_info $tmp_clc_output | grep \'Total\' | awk \'{print \$2}\'`;
			chomp $GENOMESIZE;
			printl("Predicted genome size (all reads assembly) = $GENOMESIZE bp\n");
		}
		printl("==== DONE ====\n\n");
	}else{
		printl("Could not cd to $DIR .\nExiting...\n");
		exit;	
	}
	
	#SET hash
	%dataset_hash = ();
	@{$dataset_hash{0}} = @unpaired_454_inputs;
	@{$dataset_hash{1}} = @paired_454_inputs;
	@{$dataset_hash{2}} = @unpaired_illum_inputs;
	@{$dataset_hash{3}} = @paired_illum_inputs;
	@{$dataset_hash{4}} = ();
	
	#PRE-FILTER CONTAMINATION CHECK AND SUBMISSION PREP
	if ( $presub_contam ){
		printl("\n\n==== PRESUBMISSION CONTAMINATION REMOVAL ====\n");
		#
		if ( @unpaired_454_inputs ){ 
			my %hashy = contam_rm($presub_contam, 1, "sff", 0, @unpaired_454_inputs); 
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{0}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{0}}, "$_");
					}
					
				}
			}
		}
		if ( @paired_454_inputs ){ 
			my %hashy = contam_rm($presub_contam, 1, "sff", 1, @paired_454_inputs); 
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{1}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{0}}, "$_");
					}
					
				}
			}
		}
		if ( @unpaired_illum_inputs ){ 
			my %hashy = contam_rm($presub_contam, 1, "fastq", 0, @unpaired_illum_inputs); 
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{2}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{2}}, "$_");
					}
					
				}
			}
		}
		if ( @paired_illum_inputs ){ 
			my %hashy = contam_rm($presub_contam, 1, "fastq", 1, @paired_illum_inputs); 
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{3}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{2}}, "$_");
					}
					
				}
			}
		}
		
		printl("\n\n==== READ STATS FOLLOWING PRESUBMISSION CONTAMINATION REMOVAL  ====\n");
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "002_presub_contam_rm", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("002_presub_contam_rm", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
		printl("==== DONE ====\n\n");
	}
	
	#READ CORRECTION
	if ($readCOR eq "quake"){
		printl("\n\n==== QUAKE READ CORRECTION ====\n");
		if ( (@unpaired_illum_inputs > 0) || (@paired_illum_inputs > 0) ){
			if (  @{$dataset_hash{2}} ){ @{$dataset_hash{2}} = read_correct_quake($offset, @{$dataset_hash{2}}); }
			else{ printl("No Illumina Fragments Found to Read Correct...\n"); }
			if (  @{$dataset_hash{3}} ){ @{$dataset_hash{3}} = read_correct_quake($offset, @{$dataset_hash{3}}); }
			else{ printl("No Illumina Pairs Found to Read Correct...\n"); }
		}else{
			printl"No illumina data to correct.  Moving forward...\n";
		}
		printl("==== DONE ====\n\n");
		printl("\n\n==== READ STATS FOLLOWING QUAKE READ CORRECTION  ====\n");
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "003_read_corrected", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("003_read_corrected", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
	}
	elsif($readCOR eq "allpaths"){
		printl("\n\n==== ALLPATHS READ CORRECTION (K-MER NORMALIZATION) ====\n");
		if ( (@unpaired_illum_inputs > 0) || (@paired_illum_inputs > 0) ){
			
			#FRAGMENTS
			if (  @{$dataset_hash{2}} ){ 
				@{$dataset_hash{2}} = read_correct_allpaths($offset, $ILLUMinsert, 1, @{$dataset_hash{2}});
			}
			else{ printl("No Illumina Fragments Found to Read Correct...\n"); }
			
			my $NPcount = 0;
			#PAIRS
			if (  @{$dataset_hash{3}} ){ 
				my @tmp_hash = read_correct_allpaths($offset, $ILLUMinsert, 0, @{$dataset_hash{3}});
				my @tmp_pairs = ();
				my @tmp_frags = ();
				foreach my $th (@tmp_hash){
					if (index($th, "corr.pairs") != -1){ 
						#pairs found
						printl("\nADDING TO PAIRS : $th\n");
						push (@tmp_pairs, "$th");
					}
					else{
						#fragments found
						printl("\nADDING TO FRAGMENTS : $th\n");
						push (@tmp_frags, "$th");
					}
				}
				foreach (@tmp_frags){
					push(@{$dataset_hash{2}}, "$_");
				}
				@{$dataset_hash{3}} = @tmp_pairs;
				
			}
			else{ printl("No Illumina Pairs Found to Read Correct...\n"); }
			
		}else{
			printl"No illumina data to correct.  Moving forward...\n";
		}
		printl("==== DONE ====\n\n");
		printl("\n\n==== READ STATS FOLLOWING ALLPATHS READ CORRECTION  ====\n");
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "003_read_corrected", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("003_read_corrected", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
	}
	
	#DE-DUPLICATION AND UNIQUENESS CHECK
	if ($BOOLuniqcheck == 1){
		printl("\n\n==== REMOVING EXACT DUPLICATES ====\n");
		if ( (@unpaired_illum_inputs > 0) || (@paired_illum_inputs > 0) ){
			if (  @{$dataset_hash{2}} ){ @{$dataset_hash{2}} = uniq_check_frags($offset, @{$dataset_hash{2}}); }
			else{ printl("No Illumina Fragments Found to Read Correct...\n"); }
			if (  @{$dataset_hash{3}} ){ 
				my @new_pair_illum_hash = ();
				foreach my $pairset (@{$dataset_hash{3}}){
					my @pairset_return = ();
					@pairset_return = uniq_check_pairs($offset, $pairset); 
					if ( $pairset_return[0] ){ push (@{$dataset_hash{2}}, $pairset_return[0]); }
					if ( $pairset_return[1] ){ push (@new_pair_illum_hash, $pairset_return[1]); }
				}
				@{$dataset_hash{3}} = @new_pair_illum_hash;
			}
			else{ printl("No Illumina Pairs Found to Read Correct...\n"); }
		}else{
			printl"No illumina data to correct.  Moving forward...\n";
		}
		printl("==== DONE ====\n\n");
		printl("\n\n==== READ STATS FOLLOWING UNIQUENESS CHECK  ====\n");
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "004_deduplicated", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("004_deduplicated", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
	}
	
	#DUST LOW COMPLEXITY CHECK
	if ($BOOLdustcheck == 1){
		printl("\n\n==== MASKING LOW COMPLEXITY REGIONS ====\n");
		if (  @{$dataset_hash{0}} ){ @{$dataset_hash{0}} = dust_check("sff", $offset, 0, @{$dataset_hash{0}}); }
		if (  @{$dataset_hash{1}} ){ @{$dataset_hash{1}} = dust_check("sff", $offset, $SFFinsert, @{$dataset_hash{1}}); }
		if (  @{$dataset_hash{2}} ){ @{$dataset_hash{2}} = dust_check("fastq", $offset, $SFFinsert, @{$dataset_hash{2}}); }
		if (  @{$dataset_hash{3}} ){ 
			my @tmp_hash = dust_check("fastq", $offset, $SFFinsert, @{$dataset_hash{3}}); 
				my @tmp_pairs = ();
				my @tmp_frags = ();
				foreach my $th (@tmp_hash){
					if (index($th, "unmated") != -1){ 
						#fragments found
						printl("\nADDING TO FRAGMENTS : $th\n");
						push (@tmp_frags, "$th");
					}
					else{
						#pairs found
						printl("\nADDING TO PAIRS : $th\n");
						push (@tmp_pairs, "$th");
					}
				}
				foreach (@tmp_frags){
					push(@{$dataset_hash{2}}, "$_");
				}
				@{$dataset_hash{3}} = @tmp_pairs;
		} 
		printl("\n! TOTAL BASES MASKED BY DUST = $totalDUSTchanges !\n");
		printl("==== DONE ====\n\n");
		printl("\n\n==== READ STATS FOLLOWING DUST LOW COMPLEXITY CHECK  ====\n");
		
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "005_dust_masked", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("005_dust_masked", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
	}
	
	#ILLUMINA QUALITY TRIM (default qv 20 cutoff)
	#no bool check needed
	printl("\n\n==== ILLUMINA QV TRIMMING ====\n");
	if (  @{$dataset_hash{2}} ){ @{$dataset_hash{2}} = qv_trim_frags($offset, @{$dataset_hash{2}}); }
	else{ printl("No Illumina Fragments Found to QV Trim...\n"); }
	if (  @{$dataset_hash{3}} ){ 
		my @new_pair_illum_hash = ();
		foreach my $pairset (@{$dataset_hash{3}}){
			my @pairset_return = ();
			@pairset_return = qv_trim_pairs($offset, $pairset); 
			if ( $pairset_return[0] ){ push (@{$dataset_hash{2}}, $pairset_return[0]); }
			if ( $pairset_return[1] ){ push (@new_pair_illum_hash, $pairset_return[1]); }
		}
		@{$dataset_hash{3}} = @new_pair_illum_hash;
	}
	else{ printl("No Illumina Pairs Found to QV Trim...\n"); }
	printl("==== DONE ====\n\n");
	printl("\n\n==== READ STATS FOLLOWING QV TRIM  ====\n");
	summarize_reads(%dataset_hash);
	$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "006_qv_trimmed", %dataset_hash);
	if (-s("$tmp_clc_output")){
		my @tmp_hash = ("006_qv_trimmed", "$tmp_clc_output");
		@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
		%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
		$assembly_hash_count++;
	}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
	
	#CONTAMINATION CHECK
	if ( $contam ){
		printl("\n\n==== CONTAMINATION REMOVAL ====\n");
		if ($BOOLdustcheck == 0){
			if ( @unpaired_454_inputs ){ 
				my %hashy = contam_rm($contam, 0, "sff", 0, @{$dataset_hash{0}}); 
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{0}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
						
					}
				}
			}
			if ( @paired_454_inputs ){ 
				my %hashy = contam_rm($contam, 0, "sff", 1, @{$dataset_hash{1}});  
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{1}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
						
					}
				}
			}
		}else{
			if ( @unpaired_454_inputs ){ 
				my %hashy = contam_rm($contam, 0, "fastq", 0, @{$dataset_hash{0}});  
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{0}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
						
					}
				}
			}
			if ( @paired_454_inputs ){ 
				my %hashy = contam_rm($contam, 0, "fastq", 1, @{$dataset_hash{1}});  
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{1}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
						
					}
				}
			}
		}
		if ( @unpaired_illum_inputs ){ 
			my %hashy = contam_rm($contam, 0, "fastq", 0, @{$dataset_hash{2}});  
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{2}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{2}}, "$_");
					}
					
				}
			}
		}
		if ( @paired_illum_inputs ){ 
			my %hashy = contam_rm($contam, 0, "fastq", 1, @{$dataset_hash{3}});  
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{3}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{2}}, "$_");
					}
					
				}
			}
		}
		
		printl("\n\n==== READ STATS FOLLOWING CONTAMINATION REMOVAL  ====\n");
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "007_contam_rm", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("007_contam_rm", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
		printl("==== DONE ====\n\n");
	}
	
	#ADAPTER REMOVAL
	if ( @adapters ){
		printl("\n\n==== ADAPTER REMOVAL ====\n");
		#my @return; my @return_frags;
		my %return_hash;
		if ($BOOLdustcheck == 0){
			if ( @{$dataset_hash{0}} ){ 
				my %hashy = adapter_rm($SFFinsert, "sff", 0, 1, @{$dataset_hash{0}});
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{0}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
						
					}
				}
			}
			if ( @{$dataset_hash{1}} ){ 
				my %hashy = adapter_rm($SFFinsert, "sff", 1, 1, @{$dataset_hash{1}});
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{1}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
					}
				}
			}
		}else{
			if ( @{$dataset_hash{0}} ){ 
				my %hashy = adapter_rm($SFFinsert, "fastq", 0, 2, @{$dataset_hash{0}});
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{0}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
					}
				}
			}
			if ( @{$dataset_hash{1}} ){
				my %hashy = adapter_rm($SFFinsert, "fastq", 0, 2, @{$dataset_hash{1}});
				print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
				@{$dataset_hash{1}} = @{$hashy{1}};
				if (@{$hashy{2}}){
					if (@{$hashy{2}}[0] ne "NULL"){
						foreach( @{$hashy{2}} ){
							push(@{$dataset_hash{0}}, "$_");
						}
					}
				}
			}
		}
		if ( @{$dataset_hash{2}} ){
			my %hashy = adapter_rm($SFFinsert, "fastq", 0, 3, @{$dataset_hash{2}});
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{2}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{2}}, "$_");
					}
				}
			}
		}
		if ( @{$dataset_hash{3}} ){
			my %hashy = adapter_rm($SFFinsert, "fastq", 0, 4, @{$dataset_hash{3}});
			print "FULL : @{$hashy{1}}\nPARTIAL : @{$hashy{2}}\n";
			@{$dataset_hash{3}} = @{$hashy{1}};
			if (@{$hashy{2}}){
				if (@{$hashy{2}}[0] ne "NULL"){
					foreach( @{$hashy{2}} ){
						push(@{$dataset_hash{2}}, "$_");
					}
				}
			}
		}
		printl("\n\n==== READ STATS FOLLOWING ADAPTER REMOVAL  ====\n");
		summarize_reads(%dataset_hash);
		$tmp_clc_output = assemble_clc($SFFinsert, $ILLUMinsert, "008_adapter_rm", %dataset_hash);
		if (-s("$tmp_clc_output")){
			my @tmp_hash = ("008_adapter_rm", "$tmp_clc_output");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ printl("+ ERROR : CLC assembly failed -- see logs.\nExiting...\n"); exit;}
		printl("==== DONE ====\n\n");
	}
	
	#COVERAGE FILTER
	if ($BOOLcovcheck == 1){
		printl("+ Coverage filter is not run in automaton version 1.0.  Sorry! +\n");
	}
	
	
	
	#FILTERING PIPELINE-----------------------------------------------------------------------------------------------------------STOP
	
	#...
	
	
	#12/20/12 update : output final fastq files
	my $fastfile_out = "$PREFIX" . ".post_automaton_fastq_locs.txt";
	open (LOCOUT, '>> tmp.post_automaton_fastq_locs.txt');
	print LOCOUT "\n\n! FASTQ OUTPUT FOLLOWING ALL PRE-PROCESSING:";
	print LOCOUT "\n\t! 454 FRAGMENTS: ";
	foreach(@{$dataset_hash{0}}){ print LOCOUT "\n\t\t$_"; }
	print LOCOUT "\n\t! 454 PAIRS: ";
	foreach(@{$dataset_hash{1}}){ print LOCOUT "\n\t\t$_"; }
	print LOCOUT "\n\t! ILLUMINA FRAGMENTS: ";
	foreach(@{$dataset_hash{2}}){ print LOCOUT "\n\t\t$_"; }
	print LOCOUT "\n\t! ILLUMINA PAIRS: ";
	foreach(@{$dataset_hash{3}}){ print LOCOUT "\n\t\t$_"; }
	close LOCOUT;
	runsys("mv tmp.post_automaton_fastq_locs.txt $fastfile_out");

	
	#ASSEMBLY--------------------------------------------------------------------------------------------------------------------START
		#EXPECTED COVERAGE (For Velvet)
		if ($GENOMESIZE eq "auto"){
			$GENOMESIZE = `$CLC_NGS_INSTALL/sequence_info $tmp_clc_output | grep \'Total\' | awk \'{print \$2}\'`;
			chomp $GENOMESIZE;
			printl("Predicted genome size (after all filtering) = $GENOMESIZE bp\n");
		}			
		my $bp_count = summarize_reads(%dataset_hash);
		my $exp_cov = sprintf( "%.0f", ($bp_count/$GENOMESIZE));
		printl("Predicted maximum genome coverage =  $exp_cov x\n( using Genome Span = $GENOMESIZE )\n");


		
		#VELVET
		my $velv_out;
		if ($BOOLdustcheck == 0 && (!(@adapters))){
			$velv_out = assemble_velvet($SFFinsert, $ILLUMinsert, "sff", "091_final_filtered_velvet", $exp_cov, %dataset_hash);
		}else{
			$velv_out = assemble_velvet($SFFinsert, $ILLUMinsert, "fastq", "091_final_filtered_velvet", $exp_cov, %dataset_hash);
		}
		if (-s("$velv_out")){
			my @tmp_hash = ("091_final_filtered_velvet", "$velv_out");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ 
			if ($no_assemb_switch == 1){
				printl("+ NO VELVET OUTPUT FOUND (user requested no assembly\n");
			}else{
				printl("+ ERROR : Velvet assembly failed -- see logs.\nExiting...\n"); #exit;
			}
		}
		
		#NEWBLER
		my $newb_out;
		if ($BOOLdustcheck == 0 && (!(@adapters))){
			$newb_out = assemble_newbler($SFFinsert, $ILLUMinsert, "sff", "092_final_filtered_newbler", $offset, %dataset_hash);
		}else{
			$newb_out = assemble_newbler($SFFinsert, $ILLUMinsert, "fastq", "092_final_filtered_newbler", $offset, %dataset_hash);
		}
		if (-s("$newb_out")){
			my @tmp_hash = ("092_final_filtered_newbler", "$newb_out");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ 
			if ($no_assemb_switch == 1){
				printl("+ NO NEWBLER OUTPUT FOUND (user requested no assembly\n");
			}else{
				printl("+ ERROR : Newbler assembly failed -- see logs.\nExiting...\n"); #exit;
			}
		}
		
		#assemble_celera
		my $ca_out;
		if ($BOOLdustcheck == 0 && (!(@adapters))){
			$ca_out = assemble_celera($SFFinsert, $ILLUMinsert, "sff", "093_final_filtered_celera", $offset, %dataset_hash);
		}else{
			$ca_out = assemble_celera($SFFinsert, $ILLUMinsert, "fastq", "093_final_filtered_celera", $offset, %dataset_hash);
		}
		if (-s("$ca_out")){
			my @tmp_hash = ("093_final_filtered_celera", "$ca_out");
			@{$assembly_hash{$assembly_hash_count}} = @tmp_hash;
			%{$assembly_data{$assembly_hash_count}} = %dataset_hash;
			$assembly_hash_count++;
		}else{ 
			if ($no_assemb_switch == 1){
				printl("+ NO CA OUTPUT FOUND (user requested no assembly\n");
			}else{
				printl("+ ERROR : Celera assembly failed -- see logs.\nExiting...\n"); #exit;	
			}
		}
	
#ASSEMB SWITCH		
if ($no_assemb_switch == 0){	
	
	#SUMMARIZE ASSEMBLIES--------------------------------------------------------------------------------------------------------
		chdir("$DIR") or die "Could not change directory to $DIR\n";
		
		printl("\n\nCOMPLETED ASSEMBLY LOCATIONS:\n");
		my $statcmd;
		my $newout = "./" . "$PREFIX" . "assembly_stat_cmds.txt";
		my $precount = 0;
		foreach my $z (sort keys %assembly_hash){
			printl("* $assembly_hash{$z}[0]\nLOC: $assembly_hash{$z}[1]\n");
			printl("INPUT: \n");
			printl("454 frgs    : @{$assembly_data{$z}{0}}\n");
			printl("454 pairs   : @{$assembly_data{$z}{1}}\n");
			printl("ILLUM frgs  : @{$assembly_data{$z}{2}}\n");
			printl("ILLUM pairs : @{$assembly_data{$z}{3}}\n\n");
			
			#STAT COMMANDS
			my $new_pre = "$PREFIX" . "$precount";
			my $statcmd = "$AUTOMATON_INSTALL/total_stat.pl -p $new_pre -c ";
			$precount++;
			$statcmd .= "$assembly_hash{$z}[1] -cutoff 500 -g $GENOMESIZE -k"; 
			open (PREPOUT, '>> assembly_stat_cmds.txt');
				if (  @{$assembly_data{$z}{0}} ){
					if ($BOOLdustcheck == 0 && (!(@adapters))){
						#454 input
						foreach ( @{$assembly_data{$z}{0}} ){
							$statcmd .= " -sfftitan $_";
						}
					}else{
						#illum input
						foreach ( @{$assembly_data{$z}{0}} ){
							$statcmd .= " -frg $_";
						}
					}
				}
				if (  @{$assembly_data{$z}{1}} ){
					if ($BOOLdustcheck == 0 && (!(@adapters))){
						#454 input
						foreach ( @{$assembly_data{$z}{1}} ){
							my $SFFsd = sprintf( "%.0f", ($SFFinsert/10));
							my $SFFlow = $SFFinsert - $SFFsd;
							my $SFFhigh = $SFFinsert + $SFFsd;
							$statcmd .= " -is $SFFlow $SFFhigh -sfftitan $_";
						}
					}else{
						#illum input
						foreach ( @{$assembly_data{$z}{1}} ){
							my $ILLUMsd = sprintf( "%.0f", ($ILLUMinsert/10));
							my $ILLUMlow = $ILLUMinsert - $ILLUMsd;
							my $ILLUMhigh = $ILLUMinsert + $ILLUMsd;
							$statcmd .= " -is $ILLUMlow $ILLUMhigh -ii $_";
						}
					}
				}
				if ( @{$assembly_data{$z}{2}} ){
					foreach ( @{$assembly_data{$z}{2}} ){
						$statcmd .= " -frg $_";
					}
				}
				if ( @{$assembly_data{$z}{3}} ){
					foreach ( @{$assembly_data{$z}{3}} ){
						my $ILLUMsd = sprintf( "%.0f", ($ILLUMinsert/10));
						my $ILLUMlow = $ILLUMinsert - $ILLUMsd;
						my $ILLUMhigh = $ILLUMinsert + $ILLUMsd;
						$statcmd .= " -is $ILLUMlow $ILLUMhigh -ii $_";
					}
				}
				print PREPOUT "$statcmd\n";
			close PREPOUT;
		}
		runsys("mv assembly_stat_cmds.txt $newout");
		printl("\n\nSTAT COMMANDS PRINTED TO FILE [ $newout ]\nPROCESS FINISHED SUCCESSFULLY!\n");


#ASSEMB SWITCH END
}else{
	printl("\nEND OF RUN.\nNo assembly requested.\n");
}			
	
	#END OF MAIN:
	close LOG;
	exit;
}