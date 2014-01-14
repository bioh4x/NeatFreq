#!/usr/bin/perl -w
# J. Craig Venter Institute
# NeatFreq ver 1.0
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)

use strict;
use File::Basename; #required for fileparse()
use Getopt::Std; #better command line input method
use Getopt::Long; #better command line input method

#GLOBALS
my $NEATFREQ_INSTALL = "/home/jmccorri/scripts/NF_prod_testing";
my %Opts;
my $pwd = `pwd`; chomp $pwd;
my $DIR = $pwd;
my $log=1;
my $keep = 0;
my $frag_count;
my %for_rpair_ids = ();
my %rev_rpair_ids = ();
my %HIGHprior_perbin = ();
my %ALLmate_perbin = ();
my @bin0_priority_mates = ();
my @bin0_priority_dirs = ();
my $input_status;

#OPEN LOG
if (-e("./NeatFreq.LOG.txt")){ system("rm ./NeatFreq.LOG.txt"); }
system("touch ./NeatFreq.LOG.txt");
open LOG, ">> ./NeatFreq.LOG.txt";

################################################################################################################################################################################################################################################
# printl : 
# print anything once to STDOUT and once to the LOG file -- pretty useful, eh?
################################################################################################################################################################################################################################################
sub printl {
	my $pr = shift;
	#run STDOUT print
	if ($log == 1){ print $pr; }
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
	print "USAGE :\nNeatFreq.pl -p <prefix> -r <reads.fasta> -x <cutoff> -rpairs <interleaved_pairs.fasta> -c <tallymer_counts.txt>\n\t-x defines RMKFcutoff, set to half of desired output coverage\n\t-r and -rpairs may be used alone for single library testing\n\tAll mates assumed to use inny orientation.\n";
	
	print "\t\tPrepping mates:\t\tfastq_to_fasta_qual IN.fastq OUT.fasta OUT.fasta.qual\n";
	print "\t\tgrep \'>\' E_COLI.shuffled.fasta | tr \'>\' \' \' | awk \'{print \$1}\' > allids.txt\n";
	print "\t\t$NEATFREQ_INSTALL/lib/extract_fastq_as_mates_F_R_and_frags.pl allids.txt separate\n";
	print "\t\t$NEATFREQ_INSTALL/lib/extractFasta -i OUT.fastq -idlist separate.for_mate_ids.txt -o pairs.1.fasta\n";
	print "\t\t$NEATFREQ_INSTALL/lib/extractFasta -i OUT.fastq -idlist separate.rev_mate_ids.txt -o pairs.2.fasta\n";
	print "\t\t$NEATFREQ_INSTALL/lib/extractFasta -i OUT.fastq -idlist separate.fragment_ids.txt -o fragments.fasta\n";
	
	print "\nEXAMPLE TALLYMER RUN CMDS :\n";
	print "\t$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt suffixerator -dna -pl -tis -suf -lcp -lossless -v -parts 4 -db ALL.uniqID.fna -indexname reads\n";
	print "\t$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt tallymer mkindex -mersize 19 -minocc 1 -indexname tyr-reads-minocc1 -counts -pl -esa reads\n";
	print "\t$NEATFREQ_INSTALL/lib/genometools-1.4.1/installed/bin/gt tallymer search -output qseqnum qpos counts sequence -tyr tyr-reads-minocc1 -q ALL.uniqID.fna > ALL.uniqID.19mer_counts.minocc1.txt\n";
	print "\nOPTIONAL\n\t-v [silence print to screen] (always verbose in logs)\n\t-b <bin selection method> (random or align)\n";
	print "\t\trandom = pull randomly from set (fast)\n\t\talign = run internal MSA within bins and include more low population sets (slow, -m flag required!)\n";
	print "\t-m <memory cutoff in MB> (required for use with -b random)\n";
	print "\t-z = keep all output files (bin information), normally cleaned at the end of each run\n";
	print "\t-N = provide NeatFreq install location (if not updated in script)\n\n";
	print "\nExiting.";
	
	close LOG;
	exit;
}

################################################################################################################################################################################################################################################
# logotime : 
################################################################################################################################################################################################################################################
sub logotime {
	print ".-----------------..----------------. .----------------. .----------------.  \n";
	print "| .--------------. | .--------------. | .--------------. | .--------------. |\n";  
	print "| | ____  _____  | | |  _________   | | |      __      | | |  _________   | |\n";  
	print "| ||_   \\|_   _| | | | |_   ___  |  | | |     /  \\     | | | |  _   _  |  | |\n";  
	print "| |  |   \\ | |   | | |   | |_  \\_|  | | |    / /\\ \\    | | | |_/ | | \\_|  | |\n";  
	print "| |  | |\\ \\| |   | | |   |  _|  _   | | |   / ____ \\   | | |     | |      | |\n";  
	print "| | _| |_\\   |_  | | |  _| |___/ |  | | | _/ /    \\ \\_ | | |    _| |_     | |\n";  
	print "| ||_____|\\____| | | | |_________|  | | ||____|  |____|| | |   |_____|    | |\n";  
	print "| |              | | |              | | |              | | |              | |\n";  
	print "| '--------------' | '--------------' | '--------------' | '--------------' |\n";  
	print " .----------------. .----------------. .----------------. .----------------. \n";  
	print "| .--------------. | .--------------. | .--------------. | .--------------. |\n";  
	print "| |  _________   | | |  _______     | | |  _________   | | |    ___       | |\n";  
	print "| | |_   ___  |  | | | |_   __ \\    | | | |_   ___  |  | | |  .'   '.     | |\n";  
	print "| |   | |_  \\_|  | | |   | |__) |   | | |   | |_  \\_|  | | | /  .-.  \\    | |\n";  
	print "| |   |  _|      | | |   |  __ /    | | |   |  _|  _   | | | | |   | |    | |\n";  
	print "| |  _| |_       | | |  _| |  \\ \\_  | | |  _| |___/ |  | | | \\  `-'  \\_   | |\n";  
	print "| | |_____|      | | | |____| |___| | | | |_________|  | | |  `.___.\\__|  | |\n";  
	print "| |              | | |              | | |              | | |              | |\n";  
	print "| '--------------' | '--------------' | '--------------' | '--------------' |\n";  
 	print "'----------------' '----------------' '----------------' '----------------'  \n\n";
}


################################################################################################################################################################################################################################################
# summarize_mates : 
# build a hash to retain all mate information
################################################################################################################################################################################################################################################
sub summarize_mates {
	my $rpairs = shift;
	
	#grab forward id's
	my @fpair_ids = index_reads($rpairs);
	
	#confirm size of arrays and index id's
	
	my @forward_ids;
	my @reverse_ids;
	my $z=0;
	foreach(@fpair_ids){
		if ($z==0){
			push(@forward_ids, "$_");
			$z=1;
		}else{
			push(@reverse_ids, "$_");
			$z=0;
		}
	}
	
	my $tmpsizeA = @forward_ids;
	my $tmpsizeB = @reverse_ids;
	
	if ($tmpsizeA == $tmpsizeB){
		#working -- index id's
		my $mate_count = 0;
		foreach(@forward_ids){
			my @tmparr = ($reverse_ids[$mate_count], 0); # [0=mate unused]
			$for_rpair_ids{$_} = \@tmparr;
			$mate_count++;
		}
		$mate_count = 0;
		foreach(@reverse_ids){
			my @tmparr = ($forward_ids[$mate_count], 0); # [0=mate unused]
			$rev_rpair_ids{$_} = \@tmparr;
			$mate_count++;
		}
	}
	else{
		printl("ERROR : Mates found with non-matching read ID counts:\n");
		printl("\t$tmpsizeA\t(forward)\n");
		printl("\t$tmpsizeB\t(reverse)\n\nExiting!\n");
		exit;
	}
}


################################################################################################################################################################################################################################################
# mer_counts : 
# parse the read counts and build median k-mer counts
################################################################################################################################################################################################################################################
sub mer_counts {
	my $counts_file = shift;
	
	open(GC,"< $counts_file");
	my @median_kmer_freqs;
	
	my $cur_id = 0;
	my $cur_max = 0;
	my $cur_min = 0;
	my $cur_sum;
	my @cur_array;
	my $cur_array_count = 0;
	my $zid;
	while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		my @xs = split(/\s+/, $line);
		# [0] = read id
		# [1] = coord offset
		# [2] = kmer frequency
		# [3] = sequence
		if ($xs[0] == $cur_id){
			$cur_array[$cur_array_count] = $xs[2];
			$cur_array_count++;
			if ($xs[2] > $cur_max){ $cur_max = $xs[2]; }
			if ( ($cur_min == 0) || ($xs[2] < $cur_min) ){ $cur_min = $xs[2]; }
			$cur_sum += $xs[2];
		}
		else{
			#sum
			my $arr_size = @cur_array;
			my $median_id = abs(($arr_size)/2);
			my @sorted_arr = (sort (@cur_array));
			my $median = $sorted_arr[$median_id];
			#debug 6/28
			if ($cur_array_count == 0){
				$cur_array_count = 1;
			}
			#my $mean = ($cur_sum/$cur_array_count);
			#$mean = sprintf( "%.2f", $mean );
			#printl("READ ID ( $xs[0] ) :\tMEDIAN= $median\t\|\tMEAN= $mean\t\|\tMAX= $cur_max\t\|\tMIN= $cur_min\n");
			$median_kmer_freqs[$cur_id] = $median;
			
			#iterate
			$cur_id = $xs[0];
			$cur_array_count = 0;
			$cur_max = 0;
			$cur_min = 0;
			$cur_sum = 0;
		}
	}
	#CATCH LAST SEQUENCE
	#debug 6/29
	my $arr_size = @cur_array;
	my $median_id = abs(($arr_size)/2);
	my @sorted_arr = (sort (@cur_array));
	my $median = $sorted_arr[$median_id];
	#my $mean = ($cur_sum/$cur_array_count);
	#$mean = sprintf( "%.2f", $mean );
	#printl("READ ID ( $cur_id ) :\tMEDIAN= $median\t\|\tMEAN= $mean\t\|\tMAX= $cur_max\t\|\tMIN= $cur_min\n");
	#printl("READ ID ( LAST SEQ ) :\tMEDIAN= $median\t\|\tMEAN= $mean\t\|\tMAX= $cur_max\t\|\tMIN= $cur_min\n");
	$median_kmer_freqs[$cur_id] = $median;
	#END
	
	close GC;
	return @median_kmer_freqs;	
}


################################################################################################################################################################################################################################################
# index_reads : 
# parse the sequence file and match read ids to the read counts
################################################################################################################################################################################################################################################
sub index_reads {
	my $reads_file = shift;
	
	my $arr_count = 0;
	my @read_ids = ();
	
	open(GC,"< $reads_file");
	while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		if ($line =~ m/^>/m){
			my @xs = split(/\s+/, $line);
			my $tmp = reverse($xs[0]);
			chop($tmp);
			my $tmp2 = reverse($tmp);
			$read_ids[$arr_count] = $tmp2;
			$arr_count++;
		}
	}
	close GC;
	
	return @read_ids;
}	


################################################################################################################################################################################################################################################
# index_reads_and_mates : 
# parse the sequence file and match read ids to the read counts
################################################################################################################################################################################################################################################
sub index_reads_and_mates {
	my $reads_file = shift;
	my $rpairs = shift;
	
	my $arr_count = 0;
	my @read_ids = ();
	
	#FIRST,FRAGMENTS
	if ($reads_file ne "NULL"){
		open(GC,"< $reads_file");
		while (defined(my $line = <GC>)){ #get each line
			chomp $line;
			if ($line =~ m/^>/m){
				#kill the carat
				my @xs = split(/\s+/, $line);
				my $tmp = reverse($xs[0]);
				chop($tmp);
				my $tmp2 = reverse($tmp);
				$read_ids[$arr_count] = $tmp2;
				$arr_count++;
			}
		}
		close GC;
	}
	
	my $frag_count = ($arr_count-1);
	
	#SECOND,MATES
	my @forward_ids;
	my @reverse_ids;
	my $z=0;
	open(GC,"< $rpairs");
	my $fcount=0;
	my $rcount=0;
	while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		if ($line =~ m/^>/m){
			
			#kill the carat
			my @xs = split(/\s+/, $line);
			my $tmp = reverse($xs[0]);
			chop($tmp);
			my $tmp2 = reverse($tmp);
			
			#set read_id 
			$read_ids[$arr_count] = $tmp2;
			
			#find pair ids
			if ($z==0){
				#push(@forward_ids, "$tmp2");
				$forward_ids[$fcount] = "$tmp2";
				$fcount++;
				$z=1;
			}else{
				#push(@reverse_ids, "$tmp2");
				$reverse_ids[$rcount] = "$tmp2";
				$rcount++;
				$z=0;
			}
			
			$arr_count++;
		}
	}
	close GC;
	
	#capture mate deets
	my $tmpsizeA = @forward_ids;
	my $tmpsizeB = @reverse_ids;
	
	if ($tmpsizeA == $tmpsizeB){
		#working -- index id's
		
		for (my $gg = 0; $gg < $rcount; $gg++){
			#print "$reverse_ids[$gg] $forward_ids[$gg]\n";
			my @tmparr = ("$reverse_ids[$gg]", 0); # [0=mate unused]
			@{$for_rpair_ids{$forward_ids[$gg]}} = @tmparr;
			@tmparr = ();
			@tmparr = ("$forward_ids[$gg]", 0); # [0=mate unused]
			@{$rev_rpair_ids{$reverse_ids[$gg]}} = @tmparr;
			#debug
			printl("MATE SET $gg : REV (1) $for_rpair_ids{$forward_ids[$gg]}[0] (2) $for_rpair_ids{$forward_ids[$gg]}[1] ... FOR (1) $rev_rpair_ids{$reverse_ids[$gg]}[0] (2) $rev_rpair_ids{$reverse_ids[$gg]}[1]\n");
		}
	}
	else{
		printl("ERROR : Mates found with non-matching read ID counts:\n");
		printl("\t$tmpsizeA\t(forward)\n");
		printl("\t$tmpsizeB\t(reverse)\n\nExiting!\n");
		exit;
	}
	return @read_ids;
}	


################################################################################################################################################################################################################################################
# index_reads : 
# parse the sequence file and match read ids to the read counts
################################################################################################################################################################################################################################################
sub coverage_reduce_calc{
	my $cov_in = shift;
 	my @median_mer_counts = @_;
 	
 	my $reduce;
 	my @reduce_percs;
 	
 	my $cX = 0;
 	foreach my $mmc (@median_mer_counts){
 		if ($mmc){
			my $change = $mmc - $cov_in;
			if ($change <= 0){
				$reduce = 0;
			}else{
				$reduce = ($cov_in/$mmc)*100;
				$reduce = sprintf( "%.2f", $reduce );
			}
			#printl("( $cX )\t\tF= $mmc\t(+) $change\tREDUCE%= $reduce\n");
			$reduce_percs[$cX] = $reduce;
			$cX++;
 		}
 	}
 	
 	return @reduce_percs;
	
}


################################################################################################################################################################################################################################################
# randarray
# select ids randomly from a given array
################################################################################################################################################################################################################################################
sub randarray {
	my @array = @_;
	my @rand = undef;
	my $seed = $#array + 1;
	my $randnum = int(rand($seed));
	$rand[$randnum] = shift(@array);
	while (1) {
		my $randnum = int(rand($seed));
		
		if (!(defined $rand[$randnum])){
			$rand[$randnum] = shift(@array);
		}
		last if ($#array == -1);
	}
	return @rand;
}

################################################################################################################################################################################################################################################
# rand_twoarray
# randomize 2 arrays while keeping ids the same on both
################################################################################################################################################################################################################################################
sub rand_twoarray {
	my @array = @_;
	my @rand = undef;
	my $seed = $#array + 1;
	my $randnum = int(rand($seed));
	$rand[$randnum] = shift(@array);
	while (1) {
		my $randnum = int(rand($seed));
		
		if (!(defined $rand[$randnum])){
			$rand[$randnum] = shift(@array);
		}
		last if ($#array == -1);
	}
	return @rand;
}


################################################################################################################################################################################################################################################
# prefer_mates
# select ids from array with preferential treatment of special cases
################################################################################################################################################################################################################################################
sub prefer_mates {
	my $bin = shift;
	my $ideal_perbin_reduce = shift;
	my @clustr_contents = @_;
	#%HIGHprior_perbin{bin}->@read_ids
	#%ALLmate_perbin{bin}->@read_ids
	my @output_mates = ();
	my @output_frags =();
	my @high_priority_in_clustr = ();
	my @check_last_highprior;
	my @check_last_mates;
	
	printl("-- PREFER MATES! ( bin $bin )\n");
	
	#################################
	#preclassify 1-sided mates
			printl("---- (0) classifying all low priority mates\n");
			my @maybe_mates = ();
			my @maybe_dirs = ();
			my @priority_mates = ();
			my @priority_dirs = ();
			
			#FIRST CLASSIFY MATES
			my $hpic_count = 0;
			my $skip = 0;
			my $print_onceA = 0; my $print_onceB = 0; 
			foreach(@clustr_contents){
				
				#split mates by whether or not they have already had their "friend" selected
				#9/13
				if ($skip == 0){
					#initial
					#if (exists $ALLmate_perbin{$_}){
						if ($ideal_perbin_reduce < 2){
						
							#9/13 SPECIAL CASE - minimal checks for pre-existing mates
							if (exists $for_rpair_ids{$_}){ #forward id found
								if ($for_rpair_ids{$_}[1] == 1){ #if forward id has had its reverse already selected
									#add to likely_output
									push(@priority_mates, "$_"); #first mate chosen
									push(@priority_dirs, "f");
									$skip = 1;
									printl("skipping output of more than 1 priority mate (for) for ideal reduce size $ideal_perbin_reduce\n");
								}
							}elsif(exists $rev_rpair_ids{$_}){ #reverse id found
								if ($rev_rpair_ids{$_}[1] == 1){ #if forward id has had its reverse already selected
									#add to output
										push(@priority_mates, "$_"); #first mate chosen
										push(@priority_dirs, "r");
										$skip = 1;
										printl("skipping output of more than 1 priority mate (for) for ideal reduce size $ideal_perbin_reduce\n");
								}
							}
						}else{
							#STANDARD CASE
							if (exists $ALLmate_perbin{$_}){
								if (exists $for_rpair_ids{$_}){ #forward id found
									if ($for_rpair_ids{$_}[1] == 1){ #if forward id has had its reverse already selected
										#add to likely_output
										push(@priority_mates, "$_"); #first mate chosen
										push(@priority_dirs, "f");
									}else{
										push(@maybe_mates, "$_");
										push(@maybe_dirs, "f");
									}
								}elsif(exists $rev_rpair_ids{$_}){ #reverse id found
									if ($rev_rpair_ids{$_}[1] == 1){ #if forward id has had its reverse already selected
										#add to output
											push(@priority_mates, "$_"); #first mate chosen
											push(@priority_dirs, "r");
									}else{
										push(@maybe_mates, "$_");
										push(@maybe_dirs, "r");
									}
								}
							}	
							#$skip = 2;
						}
					#}
				#9/14 update - if skip
				}else{
					printl("Skipping mate check on low ideal per bin reduce $ideal_perbin_reduce\n");	
					if ($print_onceB == 0){
						printl("Skipping high priority check on low ideal per bin reduce $ideal_perbin_reduce\n");
						$print_onceB = 1;
					}
				}
				#end 9/13
				
				#9/13
				if ($ideal_perbin_reduce >= 2){
				
					#evaluate high priority within bin
					my $tmpsize = scalar(@{$HIGHprior_perbin{$bin}});
					if ($tmpsize > 0){
						printl("highprior in clustr : ");
						for (my $z = 0; $z < $tmpsize; $z++){
							#9/7
							if (${HIGHprior_perbin{$bin}}[$z]){
								if (${HIGHprior_perbin{$bin}}[$z] eq $_){
									#push(@high_priority_in_clustr, "$_");
									$high_priority_in_clustr[$hpic_count] = "$_";
									#debug
									#printl("( $hpic_count ) $_ ");
									#printl("-- high priority found in clustr : $_\n");
									$hpic_count++;
								}
							}
						}
					}
					printl("\n");
				
				}else{
					if ($print_onceA == 0){
						printl("Skipping high priority check on low ideal per bin reduce $ideal_perbin_reduce\n");
						$print_onceA = 1;
					}
					@{$HIGHprior_perbin{$bin}} = ();
					
				}#end 9/13
				
			}
			
			my $count = scalar(@clustr_contents);
			printl("-- BIN $bin : $count (all reads in clustr)\n");
			if ($count < $ideal_perbin_reduce){
				printl("ERROR : bin $bin , count $count is less than ideal perbin reduce $ideal_perbin_reduce \nExiting.");
				exit;
			}
			$count = scalar(@{$ALLmate_perbin{$bin}});
			my $allmatecount = scalar(@priority_mates) + scalar(@maybe_mates);
			my $tmpallmatecount = $allmatecount + scalar(@high_priority_in_clustr);
			printl("-- BIN $bin : $tmpallmatecount / $count (all mates in clustr / bin)\n");
			$count = scalar(@priority_mates);
			printl("---- priority mates $count ");
			$count = scalar(@maybe_mates); 
			printl("- maybe mates $count\n");
			my $highpribincount = scalar(@{$HIGHprior_perbin{$bin}});
			my $highpricount = scalar(@high_priority_in_clustr);
			printl("-- BIN $bin : $highpricount / $highpribincount (all high priority in clustr / bin)\n");
	#################################
	
	
	#SELECT ALL IDS 
	my $id_count = 0;
	#my $highpri_skip = 0;
	while($id_count < $ideal_perbin_reduce){
		
		my $allmatecount = scalar(@priority_mates) + scalar(@maybe_mates);
		$highpricount = scalar(@high_priority_in_clustr);
		
		printl("+ ( BIN $bin ) evaluating ( $id_count / $ideal_perbin_reduce )\n");
		
		#HIGH PRIORITY FIRST
		@high_priority_in_clustr = sort @high_priority_in_clustr;
		if (($highpricount >= $ideal_perbin_reduce) && ($ideal_perbin_reduce >= 2)){
			printl("---- (1.1) all high priority export\n");
			for (my $z=0; $z < $ideal_perbin_reduce; $z+=2){
				push(@output_mates, "$high_priority_in_clustr[$z]");
				my $zz = $z+1;
				push(@output_mates, "$high_priority_in_clustr[$zz]");
				
				#debug
				printl("\tSELECTED ( $z ) : $high_priority_in_clustr[$z]\n");
				printl("\tSELECTED ( $zz ): $high_priority_in_clustr[$zz]\n");
				
				$id_count+=2;
				#update mate status
					#1
					if (exists $for_rpair_ids{$high_priority_in_clustr[$z]}){ #forward id found
						my $reverse_id = $for_rpair_ids{$high_priority_in_clustr[$z]}[0];
						$rev_rpair_ids{$reverse_id}[1] = 1;
					}elsif(exists $rev_rpair_ids{$high_priority_in_clustr[$z]}){
						my $forward_id = $rev_rpair_ids{$high_priority_in_clustr[$z]}[0];
						$for_rpair_ids{$forward_id}[1] = 1;
					}
					#2
					if (exists $for_rpair_ids{$high_priority_in_clustr[$zz]}){ #forward id found
						my $reverse_id = $for_rpair_ids{$high_priority_in_clustr[$zz]}[0];
						$rev_rpair_ids{$reverse_id}[1] = 1;
					}elsif(exists $rev_rpair_ids{$high_priority_in_clustr[$zz]}){
						my $forward_id = $rev_rpair_ids{$high_priority_in_clustr[$zz]}[0];
						$for_rpair_ids{$forward_id}[1] = 1;
					}
				#remove from array
				splice(@high_priority_in_clustr, $z, 1);
				splice(@high_priority_in_clustr, $zz, 1);
			}
			$id_count++;
			$highpricount--;
		}elsif (($highpricount > 0) && ($ideal_perbin_reduce >= 2)){
			printl("---- (1.2) all remaining high priority export\n");
			my $tmpsizeC = scalar(@high_priority_in_clustr);
			for (my $z = 0; $z < $tmpsizeC; $z=$z+2){
				push(@output_mates, "$high_priority_in_clustr[$z]");
				push(@output_mates, "$high_priority_in_clustr[($z+1)]");
				
				#debug
				#printl("\tSELECTED : $high_priority_in_clustr[$z]\n");
				#printl("\tSELECTED : $high_priority_in_clustr[($z+1)]\n");
				
				$id_count+=2;
				#update mate status
					#1
					if (exists $for_rpair_ids{$high_priority_in_clustr[$z]}){ #forward id found
						my $reverse_id = $for_rpair_ids{$high_priority_in_clustr[$z]}[0];
						$rev_rpair_ids{$reverse_id}[1] = 1;
					}elsif(exists $rev_rpair_ids{$high_priority_in_clustr[$z]}){
						my $forward_id = $rev_rpair_ids{$high_priority_in_clustr[$z]}[0];
						$for_rpair_ids{$forward_id}[1] = 1;
					}
					#2
					if (exists $for_rpair_ids{$high_priority_in_clustr[($z+1)]}){ #forward id found
						my $reverse_id = $for_rpair_ids{$high_priority_in_clustr[($z+1)]}[0];
						$rev_rpair_ids{$reverse_id}[1] = 1;
					}elsif(exists $rev_rpair_ids{$high_priority_in_clustr[($z+1)]}){
						my $forward_id = $rev_rpair_ids{$high_priority_in_clustr[($z+1)]}[0];
						$for_rpair_ids{$forward_id}[1] = 1;
					}
				#remove from array
				splice(@high_priority_in_clustr, $z, 1);
				splice(@high_priority_in_clustr, ($z+1), 1);
			}
			#empty array
			foreach(@high_priority_in_clustr){
				push(@check_last_highprior, "$_");
			}
			@high_priority_in_clustr = ();
		}
		elsif ($highpricount > 0){
			foreach(@high_priority_in_clustr){
				push(@check_last_highprior, "$_");
			}
			@high_priority_in_clustr = ();
		}
		
		#ALL MATES NEXT
		#elsif (($allmatecount > ($ideal_perbin_reduce - $id_count)) && (($ideal_perbin_reduce - $id_count) > 0)){
		elsif (($allmatecount > ($ideal_perbin_reduce - $id_count)) && (($ideal_perbin_reduce - $id_count) > 1)){
			printl("---- (2) selecting from low priority mates\n");
			
			#9/5
			my $posthi_ideal = $ideal_perbin_reduce - $id_count;
							
			#choose maximal from mates
			my $pm_size = scalar(@priority_mates);
			my $mm_size = scalar(@maybe_mates);
			my $hl_size = scalar(@check_last_highprior);
			printl("--- PMSIZE : $pm_size / MMSIZE : $mm_size\n"); #/ HPRILAST : $hl_size\n");
			
			if ($pm_size > 0){
				#if enough mates have already had their compliment taken
				if ($pm_size >= $posthi_ideal){
					printl("---- (2.1) extracting all remaining reads from previously selected pair priority subset\n");
					
					#shuffle priority_mates
					my @new_priority_mates = randarray(@priority_mates);
					my @new_priority_dirs = ();
					#update for randomized dirs
					foreach (my $aa = 0; $aa < scalar(@priority_mates); $aa++){
						foreach (my $bb = 0; $bb < scalar(@new_priority_mates); $bb++){
							if ($priority_mates[$aa] eq $new_priority_mates[$bb]){
								$new_priority_dirs[$bb] = "$priority_dirs[$aa]";
							}
						}
					}
					@priority_mates = ();
					@priority_mates = @new_priority_mates;
					@priority_dirs = ();
					@priority_dirs = @new_priority_dirs;
					
					for (my $h = 0; $h < ($posthi_ideal-1); $h++){
						#add to likely_output
							push(@output_mates, "$priority_mates[$h]");
							
							#debug
							printl("\tSELECTED : $priority_mates[$h]\n");				
							
						#update corresponding mate id status
							if ($priority_dirs[$h] eq "f"){
								my $reverse_id = $for_rpair_ids{$priority_mates[$h]}[0];
								$rev_rpair_ids{$reverse_id}[1] = 1;
							}
							elsif($priority_dirs[$h] eq "r"){
								my $forward_id = $rev_rpair_ids{$priority_mates[$h]}[0];
								$for_rpair_ids{$forward_id}[1] = 1;
							}
						#remove from pm array
							my $pm_size = scalar(@priority_mates);
							splice(@priority_mates, $h, 1);
						#update count
							$id_count++;
					}
				}
				#elsif there are still some with compliments and enough room for them all to join
				elsif($pm_size > 0){
					printl("---- (2.2) extracting all previously selected pair priority mates\n");
					for (my $h = 0; $h < $pm_size; $h++){
						#add to likely_output
							push(@output_mates, "$priority_mates[$h]");
							
							#debug
							printl("\tSELECTED : $priority_mates[$h]\n");	
							
						#update corresponding mate id status
							if ($priority_dirs[$h] eq "f"){
								my $reverse_id = $for_rpair_ids{$priority_mates[$h]}[0];
								$rev_rpair_ids{$reverse_id}[1] = 1;
							}
							elsif($priority_dirs[$h] eq "r"){
								my $forward_id = $rev_rpair_ids{$priority_mates[$h]}[0];
								$for_rpair_ids{$forward_id}[1] = 1;
							}
						#remove from pm array
							splice(@priority_mates, $h, 1);
						#update counts
							$id_count++;
					}
					#empty array
					foreach(@priority_mates){
						push(@check_last_mates, "$_");
					}
					@priority_mates = ();
				}
			}
			#elsif (($mm_size > 0){
			elsif (($mm_size > 0) && ($posthi_ideal > 1)){
				#elsif we have more than enough mates left
				if($mm_size >= $posthi_ideal){
					printl("---- (2.3) evaluating all remaining mates (compliments NOT prev. selected)\n");
					
					#shuffle maybe_mates
					#@maybe_mates = randarray(@maybe_mates);
					#shuffle maybe_mates
					my @new_maybe_mates = randarray(@maybe_mates);
					my @new_maybe_dirs = ();
					#update for randomized dirs
					foreach (my $aa = 0; $aa < scalar(@maybe_mates); $aa++){
						foreach (my $bb = 0; $bb < scalar(@new_maybe_mates); $bb++){
							if ($maybe_mates[$aa] eq $new_maybe_mates[$bb]){
								$new_maybe_dirs[$bb] = "$maybe_dirs[$aa]";
							}
						}
					}
					@maybe_mates = ();
					@maybe_mates = @new_maybe_mates;
					@maybe_dirs = ();
					@maybe_dirs = @new_maybe_dirs;
					
					
					#9/7--9/10
					for (my $h = 0; $h < ($posthi_ideal-1); $h++){
						#add to likely_output
							push(@output_mates, "$maybe_mates[$h]");
							#debug
							printl("\tSELECTED ($h) : $maybe_mates[$h]\n");	
						#update corresponding mate id status
							if ($maybe_dirs[$h] eq "f"){
								my $reverse_id = $for_rpair_ids{$maybe_mates[$h]}[0];
								$rev_rpair_ids{$reverse_id}[1] = 1;
							}
							elsif($maybe_dirs[$h] eq "r"){
								my $forward_id = $rev_rpair_ids{$maybe_mates[$h]}[0];
								$for_rpair_ids{$forward_id}[1] = 1;
							}
						#remove from mm array
							splice(@maybe_mates, $h, 1);
							$mm_size = scalar(@maybe_mates);
							printl("--- MMSIZE after reduce : $mm_size\n");
						#update counts
							$id_count++;
					}
					
				}
				#elsif we have too few mates left
				elsif($mm_size > 0){
					printl("---- (2.4) evaluating insufficient quanitity of remaining mates (compliments NOT prev. selected)\n");
					
					#shuffle maybe_mates
					#my @maybe_mates = randarray(@maybe_mates);
					
					for (my $h = 0; $h < $mm_size; $h++){
						#add to likely_output
							push(@output_mates, "$maybe_mates[$h]");
							
							#debug
							printl("\tSELECTED : $maybe_mates[$h]\n");
							
						#update corresponding mate id status
							if ($maybe_dirs[$h] eq "f"){
								my $reverse_id = $for_rpair_ids{$maybe_mates[$h]}[0];
								$rev_rpair_ids{$reverse_id}[1] = 1;
							}
							elsif ($maybe_dirs[$h] eq "r"){
								my $forward_id = $rev_rpair_ids{$maybe_mates[$h]}[0];
								$for_rpair_ids{$forward_id}[1] = 1;
							}
						#remove from mm array
							splice(@maybe_mates, $h, 1);
						#update counts
							$id_count++;
					}
					#empty array
					foreach(@maybe_mates){
						push(@check_last_mates, "$_");
					}
					@maybe_mates = ();
				}
				#move forward to fragments
				else{
					#skip to fragments
					printl("---- (2.x) waiting to proceed to fragments\n");
				}
			}else{
				#skip to fragments
				printl("---- (2.xx) waiting to proceed to fragments\n");
			}
		}
		#end all mates check
		
		#FRAGMENTS LAST
		elsif( ((scalar(@output_mates)+scalar(@output_frags)) < $ideal_perbin_reduce) && ( ($ideal_perbin_reduce-$id_count) >= 1) ){
			printl("---- (3.0) evaluating all fragments\n");
			
			#delete used mates
				# @clustr_contents vs. @output_mates
				my @diff = ();
				my %count = ();
				foreach(@clustr_contents){ $count{$_}++; }
				foreach(@output_mates){ $count{$_}++; }
				foreach my $e (keys %count){
					if ($count{$e} <= 1){ #only in clustr_contents
						push(@diff, $e);
					}
				}
				if (!(scalar(@diff) > 0)){ 
					printl("WARNING : no fragments found in comparison to mate list ( bin $bin )\n");
				}
				
				my $fr_count = scalar(@diff);
				my $hl_size = scalar(@check_last_highprior);
				printl("--- FRAGS : $fr_count / HPRILAST : $hl_size\n");
				
			#randomize fragments
				my @extractable_frgs = randarray(@diff);
			#select to fill
				my $extract_size = $ideal_perbin_reduce - (scalar(@output_mates) + scalar(@output_frags));
				#if ($extract_size <= 0){ #enough fragments exist
				#9/14
				#if ($extract_size > scalar(@extractable_frgs)){
				#9/17
				if ($extract_size < scalar(@extractable_frgs)){
					printl("---- (3.1) extracting $extract_size fragments\n");
					for (my $z=0; $z < $extract_size; $z++){
						#add to output
						push(@output_frags, "$extractable_frgs[$z]");
						
						#debug
						printl("\tSELECTED : $extractable_frgs[$z]\n");
						
						#up id count
						$id_count++;
						#no update to paired end status required!  neither of these will ever be seen again
					}
				}else{ #not enough reads!
					printl("---- (3.2) extracting all remaining fragments\n");
					printl("WARNING : not enough reads found for extraction following fragment selection ( bin $bin )\n");
					foreach(@extractable_frgs){
						#add to output
						push(@output_frags, "$_");
						#up id count
						$id_count++;
					}
					printl("---- (3.3) extracting more from high priority and single ended mates\n");
					my $lastC = 0;
					while($lastC < ($ideal_perbin_reduce-$id_count)){
						if (scalar(@check_last_highprior) > 0){
							@check_last_highprior = sort @check_last_highprior;
							push(@output_frags, "$check_last_highprior[$lastC]");
							
							#debug
							printl("\tSELECTED : $check_last_highprior[$lastC]\n");
							
							$id_count++;
							$lastC++;
						}
						my $lastD = 0;
						if (scalar(@check_last_mates) > 0){
							@check_last_mates = sort @check_last_mates;
							push(@output_frags, "$check_last_mates[$lastD]");
							
							#debug
							printl("\tSELECTED : $check_last_mates[$lastD]\n");
							
							$id_count++;
							$lastC++;
							$lastD++;
						}
						if ( (scalar(@check_last_highprior) <= 0) && (scalar(@check_last_mates) <= 0) ){
							printl("WARNING - NONE SELECTED : skipping for insufficient reads\n");
							$id_count++;
							$lastC++;
						}
					}
				}
		}
		#end fragments check
		
		else{
			#scalar(@output_mates)+scalar(@output_frags)) < $ideal_perbin_reduce) && ( ($ideal_perbin_reduce-$id_count) >1) ){
			#printl("---- (3.x) waiting to proceed to end of bin\n");
			printl("---- (3.x) end of bin selection!\n");
			printl("ideal perbin reduce : $ideal_perbin_reduce\n");
			printl("id count : $id_count\n");
			my $count1 = scalar(@output_mates); my $count2 = scalar(@output_frags);
			#printl(" scalar(@output_frags)")
			printl("output mates : $count1 / output frags : $count2\n");
			$id_count++;
		}
	}
	#END WHILE IDCOUNT LOOP
	
	#SUMMARIZE
	printl("---- (4) summarizing selection\n");
	my @all_out = ();
	foreach(@output_mates){
		push(@all_out, "$_");
	}
	foreach(@output_frags){
		push(@all_out, "$_");
	}
	
	printl("\n---- (5) exiting prefer_mates\n");
	@all_out = grep { defined } @all_out;
	return @all_out;
}


################################################################################################################################################################################################################################################
# cov_reduce
################################################################################################################################################################################################################################################
sub cov_reduce{
	my $reads_file = shift;
	my %info_hash = @_;
	
	#cleanup
	runsys ("rm ./BIN_*");
	
	#INITIALIZE BINS
	my %cov_reduce_bins;
	for (my $i=1; $i <= 99; $i++){
		@{$cov_reduce_bins{$i}} = ();
	}
	
	for my $key ( sort {$info_hash{$a}[1] <=> $info_hash{$b}[1] } keys %info_hash ){
		
		#round perc to integer for bins
		my $abs;
		if ($info_hash{$key}[2] <= 0.5 && $info_hash{$key}[2] != 0  ){ $abs = 1; }
		else{ $abs = int($info_hash{$key}[2]); }
		#print "INT = $abs\n";
		
		push(@{ $cov_reduce_bins{$abs} }, "$info_hash{$key}[0]");
	}
	
	#SUMMARIZE BINS
	printl("\n\nCOVERAGE BIN POPULATIONS:\n");
	for my $b (sort {$a<=>$b} keys %cov_reduce_bins){  
		my $blen = @{ $cov_reduce_bins{$b} };
		if ($blen > 0){ 
			my $reduced_count = int($blen * ($b/100));
			if ($reduced_count < 1){ $reduced_count = 1; }
			if ($b==0){	printl("\t100 % :\t$blen ( $blen )\n"); }#\t\t@{ $cov_reduce_bins{$b}}\n"); }
			else{ printl("\t$b % :\t$blen ( $reduced_count )\n"); }#\t\t@{ $cov_reduce_bins{$b}}\n"); }
		}
	}
	print

	#REDUCE
	printl("\n\nCOVERAGE BIN POPULATIONS:\n");
	for my $b (sort {$a<=>$b} keys %cov_reduce_bins){  
		my $blen = @{ $cov_reduce_bins{$b} };
		if ($blen > 0){ 		
			#special case
			if ($b == 0 || $b == 100){
				#file setup
				my $tmpfile = "./BIN_" . "0" . ".ids.txt";
				
				#ids
				my @tmparr = ();
				my $arr_size = @{ $cov_reduce_bins{$b} };
				if ($arr_size <= 500000){
					open ID_OUT, ">> ./BIN_0.ids.txt";
					foreach(@{ $cov_reduce_bins{$b} } ){
						push(@tmparr, "$_\n");
						print ID_OUT "$_\n"; #updated to remove carriage return - updated again 3/27
					}
					close ID_OUT;
				}else{
					my @fic_files;
					my $fic_file;
					for (my $fic = 0; $fic <= ($arr_size/500000); $fic++){
						 $fic_file = "./BIN_0_" . "$fic" . ".ids.txt";
						 open ID_OUT, "> $fic_file";
						 my $fic_count = 0;
						 for (my $x=(($fic)*500000); $x<(($fic+1)*500000); $x++) {
						 	if (${ $cov_reduce_bins{$b} }[$x]){
						 		push(@tmparr, "${ $cov_reduce_bins{$b} }[$x]");
						 		print ID_OUT "${ $cov_reduce_bins{$b} }[$x]\n"; #updated to remove carraige return
						 		$fic_count++;
						 	}
						 	#skip when it doesnt exist
						 }
						 close ID_OUT;
						 if ($fic_count > 0){
						 	push(@fic_files, "$fic_file");
						 }
					}
					printl("Summarizing Large Bin 0...\n");
					runsys("touch ./BIN_0.ids.txt");
					foreach(@fic_files){
						runsys("cat $_ >> ./BIN_0.ids.txt");
						#ERROR OUTPUT ++ RUN CMD : cat ./BIN_0_221.ids.txt; >> ./BIN_0.ids.txt
					}
				}
				
				
			}
			else{
				#extract id file
				my $reduced_count = int($blen * ($b/100));
				if ($reduced_count < 1){ $reduced_count = 1; }
				
				#file setup
				my $tmpfile = "./BIN_" . "$b" . ".ids.txt";
				
				#randomly select
				my @ran = randarray(@{ $cov_reduce_bins{$b} });
				
				#ids
				my @tmparr = ();
				open ID_OUT, ">> ./ids.tmp.txt";
				for (my $f=0; $f < $reduced_count; $f++){
					push(@tmparr, "$ran[$f]");
					print ID_OUT "$ran[$f]\n";
				}
				close ID_OUT;
				
				
				#move!
				system("mv ./ids.tmp.txt $tmpfile");
				#DEBUG 8/8/12
				if (!(-s("$tmpfile"))){
					printl("\n\nERROR: Could not update bin ( $tmpfile ) because it doesn't exist.  Permissions issue?\n\nExiting.");
					exit;
				}
			}	
		}
	}
	
	printl("Complete.\n\n");
	
	#JOIN FASTA IDS
	runsys("cat BIN_*ids.txt > ./REDUCED_COVERAGE_ID_LIST.txt");

	# DEBUG EXIT - NO MORE AS OF 1/14/14
	#printl("\n\n\nDEBUG EXIT : 3rd party software req'd\nexit.\n");
	#exit;
	
	#EXTRACT
	runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist ./REDUCED_COVERAGE_ID_LIST.txt -o ./REDUCED_COV.fasta");
}

################################################################################################################################################################################################################################################
# align_and_reduce
################################################################################################################################################################################################################################################
sub align_and_reduce{
	my $reads_file = shift;
	my $mem = shift;
	my %info_hash = @_;
	
	#cleanup
	runsys ("rm ./BIN_*");
	
	#INITIALIZE BINS
	my %cov_reduce_bins;
	for (my $i=1; $i <= 99; $i++){
		@{$cov_reduce_bins{$i}} = ();
	}
	
	for my $key ( sort {$info_hash{$a}[1] <=> $info_hash{$b}[1] } keys %info_hash ){
		#print "$key => '", join(", ", @{ $info_hash{$key} }), "'\t";
		
		#round perc to integer for bins
		my $abs;
		if ($info_hash{$key}[2] <= 0.5 && $info_hash{$key}[2] != 0  ){ $abs = 1; }
		else{ $abs = int($info_hash{$key}[2]); }
		#print "INT = $abs\n";
		
		push(@{ $cov_reduce_bins{$abs} }, "$info_hash{$key}[0]");
	}
	
	#SUMMARIZE BINS
	printl("\n\nCOVERAGE BIN POPULATIONS:\n");
	for my $b (sort {$a<=>$b} keys %cov_reduce_bins){  
		my $blen = @{ $cov_reduce_bins{$b} };
		if ($blen > 0){ 
			my $reduced_count = int($blen * ($b/100));
			if ($reduced_count < 1){ $reduced_count = 1; }
			if ($b==0){	printl("\t100 % :\t$blen ( $blen )\n\t\t@{ $cov_reduce_bins{$b}}\n"); }
			else{ printl("\t$b % :\t$blen ( $reduced_count )\n\t\t@{ $cov_reduce_bins{$b}}\n"); }
		}
	}
	
	#ALIGN AND INFLUENCE BINS
	#EXTRACT BY BIN BEFORE ALIGNMENT#####################################################	
	my @extraction_bins;
	my @init_bins;
	my @ideal_reduce_holder;
	#extract by bin
	#	my $extract_to = "BIN_" . "$retention_value" . "_" . "allreads.list.txt";
	for my $b (sort {$a<=>$b} keys %cov_reduce_bins){  
		my $blen = @{ $cov_reduce_bins{$b} }; #blen = ideal reduced count
		
			#fix 8/14
			if ($b == 0){
				$ideal_reduce_holder[$b] = $blen;
			}else{
				$ideal_reduce_holder[$b] = $blen * ($b/100);
			}
			if ($ideal_reduce_holder[$b] < 1){$ideal_reduce_holder[$b] = 1;}
		
		
		if ($blen > 0){ 		
			#special case
			if ($b == 0 || $b == 100){
				#file setup
				my $tmpfile = "./INITBIN_" . "0" . ".ids.txt";
				runsys("touch $tmpfile");
				#write to file
				open ID_OUT, ">> ./ids.tmp.txt";
				foreach my $zzz (@{ $cov_reduce_bins{$b} }){
					print ID_OUT "$zzz\n";
				}
				close ID_OUT;
				system("mv ./ids.tmp.txt $tmpfile");
				
				#EXTRACTFASTA
				my $init_tmp = "./INITBIN_" . "0";
				runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist $tmpfile -o $init_tmp");
				$init_tmp = "./INITBIN_" . "0" . "_1.fasta";
				if (!(-s("$init_tmp"))){
					printl("ERROR ENCOUNTERED : extractFasta failed before alignment within bins started\n( Bin 100, 100% retention )\n\nEXITING.\n");
				}else{
					$init_bins[0] = "$init_tmp";
				}
			}
			else{
				#extract id file
				my $reduced_count = int($blen * ($b/100));
				if ($reduced_count < 1){ $reduced_count = 1; }
				#file setup
				my $tmpfile = "./INITBIN_" . "$b" . ".ids.txt";
				runsys("touch $tmpfile");
				#write to file
				open ID_OUT, ">> ./ids.tmp.txt";
				foreach my $zzz (@{ $cov_reduce_bins{$b} }){
					print ID_OUT "$zzz\n";
				}
				close ID_OUT;
				system("mv ./ids.tmp.txt $tmpfile");
				
				#EXTRACTFASTA
				my $init_tmp = "./INITBIN_" . "$b";
				runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist $tmpfile -o $init_tmp");
				$init_tmp = "./INITBIN_" . "$b" . "_1.fasta";
				if (!(-s("$init_tmp"))){
					printl("ERROR ENCOUNTERED : extractFasta failed before alignment within bins started\n( Bin $b )\n\nEXITING.\n");
				}else{
					$init_bins[$b] = "$init_tmp";
				}
			}
		}
	}
	
	#RUN CD-HIT-EST ON BINS
	my $itr;
	my $runcmd;
	for ($itr =0; $itr <= 99; $itr++){ #FOR EACH BIN
	
		my @tmparr = ();
	
		my $outtmp = "$init_bins[$itr]" . "_cdhitest";
		if (-s("$init_bins[$itr]")){
			$runcmd = "$NEATFREQ_INSTALL/lib/cd-hit-est -i $init_bins[$itr] -o $outtmp -c 0.95 -n 10 -aS 0.4 -b 2 -G 0 -d 0 -p 1 -l 11 -M $mem";
			runsys("$runcmd");
		}else{
			printl("WARNING : Skipping cd-hit-est on non-existent ( BIN $itr )\n");
		}
		my $outtmp_out = "$outtmp" . ".clstr";

		if ((-s("$outtmp_out")) && (-s("$init_bins[$itr]"))){
			#cd-hit-est worked
			
			#BUILD HASH OF ARRAYS
			my %clstr_info = ();
			my %clstr_counts = ();
			my @clstr_components = ();
			my $num_clusters = 0;
			my $cur_id;
			open(CLSTR,"< $outtmp_out");
			while (defined(my $line = <CLSTR>)){ #get each line
				chomp $line;
				if ($line =~ m/^>/m){
					#push previous to array!
						@{$clstr_info{$cur_id}} = @clstr_components;
						my $tmpsize = @clstr_components;
						@clstr_components = ();
						$clstr_counts{$cur_id} = $tmpsize;
						$num_clusters++;
					#start new listing
						my $tmp = reverse($line); #kill the carat
						chop($tmp);
						$cur_id = reverse($tmp);  #key for %clstr_info{$cur_id}
				}
				else{
					my @xs = split(/\s+/, $line);
					my $readname = $xs[2];
					chop($readname); chop($readname); chop($readname); #kill the elipsis
					my $tmp = reverse($readname);
					chop($tmp);
					$readname = reverse($tmp);
					push(@clstr_components, "$readname");
				}
			}
			#push last entry to array!
				@{$clstr_info{$cur_id}} = @clstr_components;
				my $tmpsize = @clstr_components;
				$clstr_counts{$cur_id} = $tmpsize;
				$num_clusters++;
			close CLSTR;
			
			
			#DETERMINE CLSTR COUNT
			my $ideal_total_reduce = $ideal_reduce_holder[$itr];
			printl("(idealtotalreduce) BIN $itr : $ideal_total_reduce\n");
			my $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters;
			if ($ideal_perbin_reduce < 1){ $ideal_perbin_reduce = 1; } #make sure you pull atleast 1 from every bin
			printl("(idealperbinreduce) BIN $itr : $ideal_perbin_reduce\n");
			
			#CALCULATE IDEAL REDUCED CLSTR SIZE 
			foreach my $value (sort {$clstr_counts{$a} <=> $clstr_counts{$b}} keys %clstr_counts ){
					
				if ($clstr_counts{$value}){#check for existence
				
					#DEBUG
					printl("(clstr_counts) KEY : $value // VALUE : $clstr_counts{$value}\n");
						#printl("\t(clstr_info) VALUE : $clstr_info{$value}\n");

					#SPECIAL BIN 0 CASE
					if ( $itr == 0){
						foreach(@{$clstr_info{$value}}){
							push(@tmparr, "$_");
							#debug
							printl("\tAdding ( $_ ) to TMPARR (type 1)\n");
						}
						printl("ALLPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
					}
					elsif ( $clstr_counts{$value} <= $ideal_perbin_reduce ){
						#write all ids for extract
							foreach(@{$clstr_info{$value}}){
								push(@tmparr, "$_");
								#debug
								printl("\tAdding ( $_ ) to TMPARR (type 1)\n");
							}
							printl("RANDOMPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
						#update counts
							$ideal_total_reduce -= $clstr_counts{$value};
							$num_clusters--;
							if($num_clusters > 0){ $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters; }
							printl("(idealperbinreduce) BIN $itr : $ideal_perbin_reduce\n");		
							
					}else{
						if ($ideal_perbin_reduce >= 1){
							#write ideal_perbin_reduce num. ids for extract
								#randomly select
								my @ran = randarray(@{$clstr_info{$value}});
								#ids
								for (my $f=0; $f < $ideal_perbin_reduce; $f++){
									push(@tmparr, "$ran[$f]");
									#debug
									printl("\tAdding ( $ran[$f] ) to TMPARR (type 2)\n");
								}
								printl("RANDOMPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
							#update counts
								$ideal_total_reduce -= $ideal_perbin_reduce;
								$num_clusters--;
								if($num_clusters > 0){ $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters; }
								printl("(idealperbinreduce) BIN $itr : $ideal_perbin_reduce\n");
								
						}else{
							printl("NOPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce ( < 0 )\n");
							$num_clusters--;
							if($num_clusters > 0){ $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters; }
						}
					}	
				
				}#end new check on hash existence
			}
			
		}elsif( (!(-s("$outtmp_out"))) && (!(-s("$init_bins[$itr]"))) ){
			#EXPECTED FAIL, NO INPUT, NO BIN EXISTS
			printl("Skipping expected missing cd-hit-est info on ( BIN $itr )\n");
			next;
		}else{ 
			#FAIL ON EXISTING CDHIT OUTPUT
			printl "ERROR: cd-hit-est run failed on bin ( $itr , $init_bins[$itr] )\n";
			exit;
		}
		
		#WRITE IDS TO TEXT
			my $tmpfile = "./BIN_" . "$itr" . ".ids.txt";
			runsys("touch $tmpfile");
			if (!(-e("$tmpfile"))){
				printl("\n\nERROR: Could not create bin text file ( $tmpfile ).  Permissions issue?\n\nExiting.");
				exit;
			}
			open ID_OUT, ">> ./ids.tmp.txt";
			foreach my $zzz (@tmparr){ print ID_OUT "$zzz\n"; }
			close ID_OUT;
			#move!
			system("mv ./ids.tmp.txt $tmpfile");
			#DEBUG 8/8/12
			if (!(-s("$tmpfile"))){
				printl("\n\nERROR: Could not update bin ( $tmpfile ) because it doesn't exist.  Permissions issue?\n\nExiting.");
				exit;
			}
			
		#END FOR EACH BIN
	}
	
	printl("Complete.\n\n");
	
	#JOIN FASTA IDS
	runsys("cat BIN_*ids.txt > ./REDUCED_COVERAGE_ID_LIST.txt");
	
	#EXTRACT
	runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist ./REDUCED_COVERAGE_ID_LIST.txt -o ./REDUCED_COV.fasta");
}





################################################################################################################################################################################################################################################
# align_and_reduce_mates
################################################################################################################################################################################################################################################
sub align_and_reduce_mates{
	my $reads_file = shift;
	my $rpairs = shift;
	my $mem = shift;
	my %info_hash = @_;
	
	#cleanup
	#runsys ("mkdir BEFOREALIGN_BINS");
	#runsys ("mv ./BIN_* BEFOREALIGN_BINS/");
	
	#INITIALIZE BINS
	my %cov_reduce_bins;
	for (my $i=1; $i <= 99; $i++){
		@{$cov_reduce_bins{$i}} = ();
	}
	
	for my $key ( sort {$info_hash{$a}[1] <=> $info_hash{$b}[1] } keys %info_hash ){
		
		#round perc to integer for bins
		my $abs;
		if ($info_hash{$key}[2] <= 0.5 && $info_hash{$key}[2] != 0  ){ $abs = 1; }
		else{ $abs = int($info_hash{$key}[2]); }
		
		push(@{ $cov_reduce_bins{$abs} }, "$info_hash{$key}[0]");
	}
	
	#SUMMARIZE BINS
	printl("\n\nCOVERAGE BIN POPULATIONS:\n");
	for my $b (sort {$a<=>$b} keys %cov_reduce_bins){  
		my $blen = @{ $cov_reduce_bins{$b} };
		if ($blen > 0){ 
			my $reduced_count = int($blen * ($b/100));
			if ($reduced_count < 1){ $reduced_count = 1; }
			if ($b==0){	printl("\t100 % :\t$blen ( $blen )\n\n"); }
			else{ printl("\t$b % :\t$blen ( $reduced_count )\n"); }
		}
	}
	
	#ALIGN AND INFLUENCE BINS
	#EXTRACT BY BIN BEFORE ALIGNMENT#####################################################	
	my @extraction_bins;
	my @init_bins;
	my @ideal_reduce_holder;
	my @ALL_BIN_ALL_ID;
	my %all_mates_in_bin = ();
	for my $b (sort {$a<=>$b} keys %cov_reduce_bins){  
		
		#8/29 BIN MATES
		printl("\n---- ENTERING BIN $b ----\n");
		%all_mates_in_bin = ();
		my @HIGHpriority_mates_in_bin = ();
		my @matelog = ();
		
		my $blen = @{ $cov_reduce_bins{$b} }; #blen = ideal reduced count
		
		#fix 8/14
		if ($b == 0){
			$ideal_reduce_holder[$b] = $blen;
		}else{
			$ideal_reduce_holder[$b] = $blen * ($b/100);
		}
		if ($ideal_reduce_holder[$b] < 1){$ideal_reduce_holder[$b] = 1;}
		
		
		if ($blen > 0){ 		
			#special case
			if ($b == 0 || $b == 100){
				#file setup
				my $tmpfile = "./INITBIN_" . "0" . ".ids.txt";
				runsys("touch $tmpfile");
				#write to file
				open ID_OUT, ">> ./ids.tmp.txt";
				foreach my $zzz (@{ $cov_reduce_bins{$b} }){
					print ID_OUT "$zzz\n";
				}
				close ID_OUT;
				system("mv ./ids.tmp.txt $tmpfile");
				
				
				my $init_frg = "./INITBINfrg_" . "0";
				my $init_prs = "./INITBINprs_" . "0";
				if ($reads_file ne "NULL"){
					runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist $tmpfile -o $init_frg");
				}
				runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $rpairs -idlist $tmpfile -o $init_prs");
				$init_frg = "./INITBINfrg_" . "0" . "_1.fasta";
				$init_prs = "./INITBINprs_" . "0" . "_1.fasta";
				my $init_tmp = "./INITBIN_" . "0" . ".fasta";
				if (-s "$init_frg"){ runsys("cat $init_frg >> $init_tmp"); }
				else{ printl("WARNING : only pairs found ( bin 0 )\n"); }
				if (-s "$init_prs"){ runsys("cat $init_prs >> $init_tmp"); }
				else{ printl("WARNING : only fragments found ( bin 0 )\n"); }
				if (!(-s("$init_tmp"))){
					printl("WARNING : extractFasta failed before alignment within bins started\n( Bin 0 )\n\nEXITING.\n");
				}else{
					$init_bins[$b] = "$init_tmp";
				}
				
			}
			else{
				#extract id file
				my $reduced_count = int($blen * ($b/100));
				if ($reduced_count < 1){ $reduced_count = 1; }
				#file setup
				my $tmpfile = "./INITBIN_" . "$b" . ".ids.txt";
				runsys("touch $tmpfile");
				#write to file
				open ID_OUT, ">> ./ids.tmp.txt";
				
				#ITERATE THROUGH BINS
				foreach my $zzz (@{ $cov_reduce_bins{$b} }){
					print ID_OUT "$zzz\n";
					
					#frst check for all mates
					if (exists $for_rpair_ids{$zzz}){
						$all_mates_in_bin{$zzz} = "f";
						push(@matelog, "$zzz");
					}elsif(exists $rev_rpair_ids{$zzz}){
						$all_mates_in_bin{$zzz} = "r";
						push(@matelog, "$zzz");
					}
					
				}
				#END ITERATE THROUGH BINS
				close ID_OUT;
				system("mv ./ids.tmp.txt $tmpfile");
				
				#8/30 -- COMPARE TO MATE PAIR LISTING ############
				#second check for for/rev matches
				my $tmpcount = keys %all_mates_in_bin;
				if ($tmpcount > 0){
					printl("---- mates from bin $b : $tmpcount\n");
					
					#compare all
					foreach my $value (sort(keys %all_mates_in_bin)){
						if ($all_mates_in_bin{$value}){#check for existence
							if ((exists $all_mates_in_bin{$value}) && ($all_mates_in_bin{$value} eq "f")){
								if (exists $for_rpair_ids{$value}){
									if (exists($for_rpair_ids{$value}[0])){
										
										#9/4updatecheckwithinbin
										if (exists($all_mates_in_bin{$for_rpair_ids{$value}[0]})){
										
											push (@HIGHpriority_mates_in_bin, "$value");
											push (@HIGHpriority_mates_in_bin, "$for_rpair_ids{$value}[0]");
											#DEBUG
											#printl("\tadding to HIGHPRI mates (found) : $value\n");
											#printl("\tadding to HIGHPRI mates (compliment) : $for_rpair_ids{$value}[0]\n");
									
										#9/4_end_updatecheckwithinbin
										}
									
									}else{
										printl("ERROR : coresponding mate could not be located for forward mate id $value\nExiting.\n\n");
										exit;
									}
								}else{
									printl("ERROR : no forward mate info found for id $value\nExiting.\n\n");
									exit;
								}
							}elsif ((exists $all_mates_in_bin{$value}) && ($all_mates_in_bin{$value} eq "r")){
								if (exists $rev_rpair_ids{$value}){
									if (exists($rev_rpair_ids{$value}[0])){
										
										#9/4updatecheckwithinbin
										if (exists($all_mates_in_bin{$rev_rpair_ids{$value}[0]})){
										
											push (@HIGHpriority_mates_in_bin, "$value");
											push (@HIGHpriority_mates_in_bin, "$rev_rpair_ids{$value}[0]");
											#DEBUG
											#printl("\tadding to HIGHPRI mates (found) : $value\n");
											#printl("\tadding to HIGHPRI mates (compliment) : $rev_rpair_ids{$value}[0]\n");
										
										#9/4_end_updatecheckwithinbin
										}
										
									}else{
										printl("ERROR : coresponding mate could not be located for reverse mate id $value\nExiting.\n\n");
										exit;
									}
								}else{
									printl("ERROR : no reverse mate info found for id $value\nExiting.\n\n");
									exit;
								}
							}
						}else{
							if (exists $all_mates_in_bin{$value}){
								printl("ERROR : no direction found in all_mates_in_bin arr $value\nExiting.\n\n");
								exit;
							}
							else{
								printl("ERROR : id $value does not appear in all_mates_in_bin!\nExiting.\n\n");
								exit;
							}
						}
					}#end for
					#remove highpriority dupes
					printl("Deduplicating High Priorirty Mates of Bin $b ...\n");
					$tmpcount = @HIGHpriority_mates_in_bin;
					#printl("---- size before dedupe from bin $b : $tmpcount\n");
					my %unique = ();
					foreach my $item (@HIGHpriority_mates_in_bin)
					{
						$unique{$item} ++;
					}
					my @HIGHpriority_mates_in_bin = keys %unique;
					$tmpcount = @HIGHpriority_mates_in_bin;
					printl("---- high priority mates from bin $b, complete F/R : $tmpcount\n");
					#ADD TO MATE TRACKING
					#9/6
					@{$HIGHprior_perbin{$b}} = @HIGHpriority_mates_in_bin;
					@{$ALLmate_perbin{$b}} = @matelog;
				}else{
					printl("---- mates from bin $b : 0\n");
				}
				##################################################
				

				my $init_frg = "./INITBINfrg_" . "$b";
				my $init_prs = "./INITBINprs_" . "$b";
				if ($reads_file ne "NULL"){
					runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist $tmpfile -o $init_frg");
				}
				runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $rpairs -idlist $tmpfile -o $init_prs");
				$init_frg = "./INITBINfrg_" . "$b" . "_1.fasta";
				$init_prs = "./INITBINprs_" . "$b" . "_1.fasta";
				my $init_tmp = "./INITBIN_" . "$b" . ".fasta";
				if (-s "$init_frg"){ runsys("cat $init_frg >> $init_tmp"); }
				else{ printl("WARNING : only pairs found ( bin $b )\n"); }
				if (-s "$init_prs"){ runsys("cat $init_prs >> $init_tmp"); }
				else{ printl("WARNING : only fragments found ( bin $b )\n"); }
				if (!(-s("$init_tmp"))){
					printl("WARNING : extractFasta failed before alignment within bins started\n( Bin $b )\n\nEXITING.\n");
				}else{
					$init_bins[$b] = "$init_tmp";
				}
				
				#CLEANUP
				if ($keep == 0){
					runsys("rm -r $init_frg");
					runsys("rm -r $init_prs");
				}
			}
		}
	}
	#FINISH ITERATING THROUGH BINS, (EXTRACTING TO INITBINS AND COUNTING PAIRS)
	
	#RUN CD-HIT-EST ON BINS
	my $itr;
	my $runcmd;
	for ($itr =0; $itr <= 99; $itr++){ #FOR EACH BIN
		
		my @tmparr = ();
		
		my $outtmp = "$init_bins[$itr]" . "_cdhitest";
		if (-s("$init_bins[$itr]")){
			$runcmd = "/$NEATFREQ_INSTALL/lib/cd-hit-est -i $init_bins[$itr] -o $outtmp -c 0.95 -n 10 -aS 0.4 -b 2 -G 0 -d 0 -p 1 -l 11 -M $mem";
			runsys("$runcmd");
		}else{
			printl("WARNING : Skipping cd-hit-est on non-existent ( BIN $itr )\n");
		}
		my $outtmp_out = "$outtmp" . ".clstr";
		
		if ((-s("$outtmp_out")) && (-s("$init_bins[$itr]"))){
			#cd-hit-est worked
			
			#BUILD HASH OF ARRAYS
			my %clstr_info = ();
			my %clstr_counts = ();
			my @clstr_components = ();
			my $num_clusters = 0;
			my $cur_id = "NULL";
			open(CLSTR,"< $outtmp_out");
			while (defined(my $line = <CLSTR>)){ #get each line
				chomp $line;
				#9/6
				if ((($line =~ m/^>/m)) && ($cur_id eq "NULL")){
					#start new listing
					my $tmp = reverse($line); #kill the carat
					chop($tmp);
					$cur_id = reverse($tmp);  #key for %clstr_info{$cur_id}
					#debug
					#printl("new cur id $cur_id\n");
				}
				elsif ($line =~ m/^>/m){
					#push previous to array!
					#9/6
					my $tmpsize = scalar(@clstr_components);
					#printl("  original size of finished cluster : $tmpsize\n");
					if ($tmpsize > 1){
						#check for dupes caused by lines with asterisks in cdhitest output
						#clean up dupes
						my @diff = ();
						my %count = ();
						foreach(@clstr_components){ $count{$_}++; }
						foreach my $e (keys %count){
							if ($count{$e} <= 1){ #only in clustr_contents
								push(@diff, $e);
							}
						}
						@{$clstr_info{$cur_id}} = @diff;
						$tmpsize = scalar(@diff);
					}else{
						@{$clstr_info{$cur_id}} = @clstr_components;
					}
					
					#$clstr_info{$cur_id} = \@clstr_components;
					@clstr_components = ();
					$clstr_counts{$cur_id} = $tmpsize;
					#debug
					#printl("  size of finished cluster : $tmpsize\n");
					$num_clusters++;
					#start new listing
					my $tmp = reverse($line); #kill the carat
					chop($tmp);
					$cur_id = reverse($tmp);  #key for %clstr_info{$cur_id}
					#debug
					#printl("new cur id $cur_id\n");
				}
				else{
					my @xs = split(/\s+/, $line);
					my $readname = $xs[2];
					chop($readname); chop($readname); chop($readname); #kill the elipsis
					my $tmp = reverse($readname);
					chop($tmp);
					$readname = reverse($tmp);
					push(@clstr_components, "$readname");
				}
			}
			#push last entry to array!
			my $tmpsize = scalar(@clstr_components);
			if ($tmpsize > 1){
				#check for dupes caused by lines with asterisks in cdhitest output
				#clean up dupes
				my @diff = ();
				my %count = ();
				foreach(@clstr_components){ $count{$_}++; }
				foreach my $e (keys %count){
					if ($count{$e} <= 1){ #only in clustr_contents
						push(@diff, $e);
					}
				}
				#9/20x
				@diff = grep { defined } @diff;
				
				@{$clstr_info{$cur_id}} = @diff;
				$tmpsize = scalar(@diff);
			}else{
				#9/20x
				@clstr_components = grep { defined } @clstr_components;
				@{$clstr_info{$cur_id}} = @clstr_components;
			}
			#debug
			#printl("  size of finished cluster : $tmpsize\n");
			$clstr_counts{$cur_id} = $tmpsize;
			$num_clusters++;
			close CLSTR;
			
			
			#DETERMINE CLSTR COUNT
			my $ideal_total_reduce = $ideal_reduce_holder[$itr];
			printl("\n\n(idealtotalreduce) BIN $itr : $ideal_total_reduce\n");
			#my $ideal_perbin_reduce = abs($ideal_total_reduce/$num_clusters);
			my $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters;
			if ($ideal_perbin_reduce < 1){ $ideal_perbin_reduce = 1; } #make sure you pull atleast 1 from every bin
			printl("(idealperbinreduce) BIN $itr : $ideal_perbin_reduce\n");
			
			#CALCULATE IDEAL REDUCED CLSTR SIZE 
			#foreach my $value (sort {$clstr_counts{$a} cmp $clstr_counts{$b}} keys %clstr_counts ){
			foreach my $value (sort {$clstr_counts{$a} <=> $clstr_counts{$b}} keys %clstr_counts ){
				
				if (exists $clstr_counts{$value}){#check for existence
					
					#DEBUG
					printl("-- (clstr_counts) KEY : $value // VALUE : $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
					#9/6
					#printl("\t(clstr_info) VALUE : @{$clstr_info{$value}}\n");
					
					#SPECIAL BIN 0 CASE
					if ( $itr == 0){
						foreach(@{$clstr_info{$value}}){
							push(@tmparr, "$_");
							#debug
							printl("\t-- Adding ( $_ ) to TMPARR (type 0)\n");
							
							
							#update 9/5 for bin 0 1-sided mates
							
							#update 9/13 for no false checks
							#if (exists $all_mates_in_bin{$_}){
							
								#update 9/18 fix correct mate ids
								if (exists $for_rpair_ids{$_}){ #forward id found
									my $reverse_id = $for_rpair_ids{$_}[0];
									$rev_rpair_ids{$reverse_id}[1] = 1;
									
									#DEBUG
									printl("FOR FOUND/REV SET : FOR (1) $for_rpair_ids{$_}[0] (2) $for_rpair_ids{$_}[1] ... REV (1) $rev_rpair_ids{$reverse_id}[0] (2) $rev_rpair_ids{$reverse_id}[1]\n");
									
								}elsif(exists $rev_rpair_ids{$_}){ #reverse id found
								 	my $forward_id = $rev_rpair_ids{$_}[0];
									$for_rpair_ids{$forward_id}[1] = 1;
									
									#DEBUG
									printl("REV FOUND/FOR SET : FOR (1) $for_rpair_ids{$forward_id}[0] (2) $for_rpair_ids{$forward_id}[1] ... REV (1) $rev_rpair_ids{$_}[0] (2) $rev_rpair_ids{$_}[1]\n");
									
								}else{
									printl("WARNING : no for or rev id match found\n");
								}
								
							#}
							#end 9/13
							
						}
						printl("\tALLPULLBINZERO ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
					}
					#BIN >= 1
					elsif ( $clstr_counts{$value} <= $ideal_perbin_reduce ){
						#write all ids for extract
						foreach(@{$clstr_info{$value}}){
							push(@tmparr, "$_");
							printl("\t-- Adding ( $_ ) to TMPARR (type 1)\n");
							
							#update 9/5 for tracking 1-sided mates
							
							#update 9/13 for no false checks
							#if (exists $all_mates_in_bin{$_}){
								
								#update 9/18 fix correct mate ids
								if (exists $for_rpair_ids{$_}){ #forward id found
									my $reverse_id = $for_rpair_ids{$_}[0];
									$rev_rpair_ids{$reverse_id}[1] = 1;
									
									#DEBUG
									printl("FOR FOUND/REV SET : FOR (1) $for_rpair_ids{$_}[0] (2) $for_rpair_ids{$_}[1] ... REV (1) $rev_rpair_ids{$reverse_id}[0] (2) $rev_rpair_ids{$reverse_id}[1]\n");
									
								}elsif(exists $rev_rpair_ids{$_}){ #reverse id found
								 	my $forward_id = $rev_rpair_ids{$_}[0];
									$for_rpair_ids{$forward_id}[1] = 1;
									
									#DEBUG
									printl("REV FOUND/FOR SET : FOR (1) $for_rpair_ids{$forward_id}[0] (2) $for_rpair_ids{$forward_id}[1] ... REV (1) $rev_rpair_ids{$_}[0] (2) $rev_rpair_ids{$_}[1]\n");
									
								}else{
									printl("WARNING : no for or rev id match found\n");
								}
							
							#}
							#end 9/13
							
						}
						printl("\tALLPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
						#update counts
						$ideal_total_reduce -= $clstr_counts{$value};
						$num_clusters--;
						if($num_clusters > 0){ $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters; }
						printl("(idealperbinreduce) BIN $itr : $ideal_perbin_reduce\n");		
						
					}else{
						if ($ideal_perbin_reduce >= 1){
							printl("\t-- DYNAMICMATEPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce\n");
							
							#PREFER_MATES!
							my @ran = prefer_mates($itr, $ideal_perbin_reduce, @{$clstr_info{$value}});

							foreach(@ran){
								push(@tmparr, "$_");
								printl("\t-- Adding ( $_ ) to TMPARR (type 2)\n");
								
								#update 9/13 for no false checks
								#if (exists $all_mates_in_bin{$_}){
									
									#update 9/18 fix correct mate ids
									if (exists $for_rpair_ids{$_}){ #forward id found
										my $reverse_id = $for_rpair_ids{$_}[0];
										$rev_rpair_ids{$reverse_id}[1] = 1;
										
										#DEBUG
										printl("FOR FOUND/REV SET : FOR (1) $for_rpair_ids{$_}[0] (2) $for_rpair_ids{$_}[1] ... REV (1) $rev_rpair_ids{$reverse_id}[0] (2) $rev_rpair_ids{$reverse_id}[1]\n");
										
									}elsif(exists $rev_rpair_ids{$_}){ #reverse id found
									 	my $forward_id = $rev_rpair_ids{$_}[0];
										$for_rpair_ids{$forward_id}[1] = 1;
										
										#DEBUG
										printl("REV FOUND/FOR SET : FOR (1) $for_rpair_ids{$forward_id}[0] (2) $for_rpair_ids{$forward_id}[1] ... REV (1) $rev_rpair_ids{$_}[0] (2) $rev_rpair_ids{$_}[1]\n");
										
									}else{
										printl("WARNING : no for or rev id match found\n");
									}
								
								#}
								#end 9/13
								
							}
							#update counts
							$ideal_total_reduce -= $ideal_perbin_reduce;
							$num_clusters--;
							if($num_clusters > 0){ $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters; }
							printl("(idealperbinreduce) BIN $itr : $ideal_perbin_reduce\n");
							
						}else{
							printl("\t-- NOPULL ( BIN $itr ) :  total= $clstr_counts{$value} // ideal= $ideal_perbin_reduce ( < 0 )\n\n");
							$num_clusters--;
							if($num_clusters > 0){ $ideal_perbin_reduce = $ideal_total_reduce/$num_clusters; }
						}
					}	
					
				}#end new check on hash existence
			}
			
		}elsif( (!(-s("$outtmp_out"))) && (!(-s("$init_bins[$itr]"))) ){
			#EXPECTED FAIL, NO INPUT, NO BIN EXISTS
			printl("Skipping expected missing cd-hit-est info on ( BIN $itr )\n");
			next;
		}else{ 
			#FAIL ON EXISTING CDHIT OUTPUT
			printl "ERROR: cd-hit-est run failed on bin ( $itr , $init_bins[$itr] )\n";
			exit;
		}
		
		#WRITE IDS TO TEXT
		my $tmpfile = "./BIN_" . "$itr" . ".ids.txt";
		runsys("touch $tmpfile");
		if (!(-e("$tmpfile"))){
			printl("\n\nERROR: Could not create bin text file ( $tmpfile ).  Permissions issue?\n\nExiting.");
			exit;
		}
		open ID_OUT, ">> ./ids.tmp.txt";
		foreach my $zzz (@tmparr){ 
			print ID_OUT "$zzz\n"; 
			push(@ALL_BIN_ALL_ID, "$zzz");
		}
		close ID_OUT;
		#move!
		system("mv ./ids.tmp.txt $tmpfile");
		#DEBUG 8/8/12
		if (!(-s("$tmpfile"))){
			printl("\n\nERROR: Could not update bin ( $tmpfile ) because it doesn't exist.  Permissions issue?\n\nExiting.");
			exit;
		}
		
		#END FOR EACH BIN
	}
	
	printl("Complete.\n\n");
	
	#JOIN FASTA IDS
	runsys("cat BIN_*ids.txt > ./REDUCED_COVERAGE_ID_LIST.txt");
	
	
	
	#PREP FOR EXTRACT
	@ALL_BIN_ALL_ID = sort @ALL_BIN_ALL_ID;
	if (scalar(@ALL_BIN_ALL_ID) > 0){
		printl("-- printing out pairs\n");
		my $fic_file = "./ALL.pairs.ids.txt";
		my $ficfrg_file = "./ALL.frags.ids.txt";
		system("touch $fic_file");
		open ID_OUT, ">> $fic_file";
		open IDfrg_OUT, ">> $ficfrg_file";
		my $mates_as_mates = 0;
		my $mates_as_frags = 0;
		my $frags_as_frags = 0;
		
		foreach(@ALL_BIN_ALL_ID){
			
			#if (exists $all_mates_in_bin{$_}){
				
				#check mate status
				if (exists $for_rpair_ids{$_}){
					if ($for_rpair_ids{$_}[1] == 1){ 
						#forward id with complimentary mate
						print ID_OUT "$_\n";
						$mates_as_mates++;
					}else{
						print IDfrg_OUT "$_\n";
						$mates_as_frags++;
					}
				}elsif(exists $rev_rpair_ids{$_}){
					if ($rev_rpair_ids{$_}[1] == 1){
						#reverse id with complimentary mate
						print ID_OUT "$_\n";
						$mates_as_mates++;
					}else{
						print IDfrg_OUT "$_\n";
						$mates_as_frags++;
					}
				}else{
					print IDfrg_OUT "$_\n";
					$frags_as_frags++;
				}
			
			#}
		}
		close ID_OUT;
		close IDfrg_OUT;
		printl("-- Mates added as mates:          $mates_as_mates\n");
		printl("-- Mates added as fragments:      $mates_as_frags\n");
		printl("-- Fragments added as fragments:  $frags_as_frags\n");
	}
	
	#EXTRACT
	#runsys("/usr/local/devel/BCIS/assembly/tools/extractFasta -i $reads_file -idlist ./REDUCED_COVERAGE_ID_LIST.txt -o ./REDUCED_COV.fragments.fasta");
	if ($reads_file ne "NULL"){
		runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $reads_file -idlist ./ALL.frags.ids.txt -o ./REDUCED_COV.fragments.fasta");
	}
	#9/20
	#runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $rpairs -idlist ./ALL.pairs.ids.txt -o ./REDUCED_COV.pairs_as_fragments.fasta");
	runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $rpairs -idlist ./ALL.frags.ids.txt -o ./REDUCED_COV.pairs_as_fragments.fasta");
	#runsys("/usr/local/devel/BCIS/assembly/tools/extractFasta -i $rpairs -idlist ./REDUCED_COVERAGE_ID_LIST.txt -o ./REDUCED_COV.pairs.fasta");
	runsys("$NEATFREQ_INSTALL/lib/extractFasta -i $rpairs -idlist ./ALL.pairs.ids.txt -o ./REDUCED_COV.pairs.fasta");
}






################################################################################################################################################################################################################################################
# MAIN : 
################################################################################################################################################################################################################################################
MAIN : {	
	my ($reads_file, $counts_file, $prefix, $cov_in, $bin_extract_type, $mem, $rpairs, $NEW_NEATFREQ_INSTALL);
	my $pairstatus = 0;
	my $status = GetOptions(\%Opts, "help!", "h!", "v!", "z!", 'b=s'=> \$bin_extract_type, 'N=s'=> \$NEW_NEATFREQ_INSTALL,'m=s'=> \$mem, 'r=s'=> \$reads_file, 'c=s'=> \$counts_file, 'p=s'=> \$prefix, 'x=s'=> \$cov_in, 'rpairs=s'=> \$rpairs);
	
	#PARSE INPUT------------------------------------------------------------------------------------------------------------------START
	print "\nNeatFreq ver 1.0\n";
	
	if ( exists $Opts{help}){ printl("Help requested:\n"); outhelp(); }	
	if ( exists $Opts{h}){ printl("Help requested:\n"); outhelp(); }
	if ( exists $Opts{v}){ 
		printl("Print of logs to screen disabled.\n");
		$log = 0;
	}else{
		logotime();
	}
	if ( exists $Opts{z}){ 
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
	if ( $prefix ){
		printl("Using User-Input Prefix : $prefix\n");
	}else{
		printl("Using Default Prefix : OUT\n");
		$prefix = "OUT";
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
	
	if (-s($counts_file)){
		printl("Using input counts file : $counts_file\n");
	}else{
		printl("No read file input or counts file is 0 bytes in size. ( $counts_file ) EXITING!\n\n");
		outhelp();
	}
	
	if ($reads_file && (-s("$reads_file")) ){
		printl("Using input reads file : $reads_file\n");
	}else{
		if ($reads_file){ printl("ERROR:\tFragment-only read file input or read file is 0 bytes in size. ( $reads_file )\nEXITING!\n"); exit; }
		else { printl("WARNING:\tNo fragment-only read file input ( $reads_file )!\n"); }
		# outhelp();
	}
	
	if ( $rpairs && (-s("$rpairs"))){
		printl("Paired input reads input: $rpairs\n");	
		$pairstatus=1;
	}else{
		printl("WARNING:\tNo paired read file input or reads files are 0 bytes in size. ( $rpairs !\n");
		# outhelp();
	}
	
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
		outhelp();
	}

	#INDEX READS
	my @read_ids;
	if ($pairstatus == 0){
		@read_ids = index_reads($reads_file);
	}
	elsif ($pairstatus == 1){
		if ($input_status == 2){
			@read_ids = index_reads_and_mates($reads_file, $rpairs);
		}elsif ($input_status == 1){
			@read_ids = index_reads_and_mates("NULL", $rpairs);
		#}elsif ($input_status == 0){
		#	@read_ids = index_reads_and_mates($reads_file, $rpairs);
		#}
		}else{
			printl("Error encountered during parse of input reads.\n\nExiting...\n");
			exit;
		}
		#summarize_mates($rpairs);
	}
	
	#ITERATE THROUGH COUNTS AND BUILD MEDIAN KMER CALCULATIONS
	printl("\nBuilding mer counts...\n");
	my @median_mer_counts = mer_counts($counts_file);
	
	#INITAL SUMMARY
	my @reduce_percs = coverage_reduce_calc($cov_in, @median_mer_counts);
	
	#MERGE INFO
	printl("\nBuilding combined hash...\n");
	my %info_hash;
	my %no_reduce;
	#FRAGMENTS
	my $len = @read_ids;
	for (my $z = 0; $z < $len; $z++){
		if ( ($read_ids[$z]) && ($median_mer_counts[$z]) && (($reduce_percs[$z]) || $reduce_percs[$z] == 0) ){
			my @tmp_arr = ($read_ids[$z], $median_mer_counts[$z], $reduce_percs[$z]);
			$info_hash{$z} = \@tmp_arr;
		}else{
			printl("\nError : data missing\n");
			printl("\t$read_ids[$z]\n\t$median_mer_counts[$z]\n\t$reduce_percs[$z]\nEXITING!\n\n");
			exit;
		}
	}
	printl("Done.\n\n");
	
	
	#REDUCE!
	if ($bin_extract_type eq "random"){
		if ($pairstatus==0){
			cov_reduce($reads_file, %info_hash);
		}else{
			printl("Merging reads for downstream random analysis...\n");
			my $new_all_reads = "./" . "$prefix" . "_allreads.fasta";
			if ($input_status==2){
				runsys("cat $reads_file > $new_all_reads");
				runsys("cat $rpairs >> $new_all_reads");
			}
			elsif ($input_status==1){
				runsys("cp $rpairs $new_all_reads");
			}
			else{ #INPUT_STATUS=0
				runsys("cp $reads_file $new_all_reads");
			}
			#runsys("cat $reads_file > $new_all_reads");
			#runsys("cat $rpairs >> $new_all_reads");
			$reads_file = $new_all_reads;
			#if(-s("$reads_file")){
			#	cov_reduce($reads_file, %info_hash);
			if(-s("$new_all_reads")){
				cov_reduce($new_all_reads, %info_hash);
			}else{ 
				printl("ERROR : Failed to merge fragments and mates to one file for random selection.\n\nExiting.\n"); 
				exit; 
			}
			$pairstatus = 0;
		}
	}
	elsif($bin_extract_type eq "align"){
		if ($pairstatus==0){
			align_and_reduce($reads_file, $mem, %info_hash);
		}else{
			if ($input_status == 2){
				align_and_reduce_mates($reads_file, $rpairs, $mem, %info_hash);
			}elsif ($input_status == 1){
				align_and_reduce_mates("NULL", $rpairs, $mem, %info_hash);
			}
			
		}
	}
	
	#CLEANUP
	runsys("mkdir KMER_BINS KMER_BINS_mp_frg");
	runsys("mv ./BIN_* KMER_BINS/");
	runsys("mv ./BINfrg_* KMER_BINS_mp_frg/");
	runsys("mv ./BINprs_* KMER_BINS_mp_frg/");
	
	#DEBUG
	#exit;
	
	my $outfile = "./" . "$prefix" . "reducedCOV_" . "$cov_in" . "X" . ".fasta";
	my $outfile_frg = "./" . "$prefix" . "reducedCOV_" . "$cov_in" . "X.fragments.fasta";
	my $outfile_prs = "./" . "$prefix" . "reducedCOV_" . "$cov_in" . "X.pairs.fasta";
	if (-s("./REDUCED_COV.fasta_1.fasta")){
		runsys("mv ./REDUCED_COV.fasta_1.fasta $outfile");
	}else{
		runsys("touch $outfile_frg");
		if (-s("./REDUCED_COV.fragments.fasta_1.fasta")){ 
			runsys("cat ./REDUCED_COV.fragments.fasta_1.fasta > $outfile_frg"); 
			runsys("cat ./REDUCED_COV.pairs_as_fragments.fasta_1.fasta >> $outfile_frg");
		}else{
			runsys("mv ./REDUCED_COV.pairs_as_fragments.fasta_1.fasta $outfile_frg");
		}
		runsys("mv ./REDUCED_COV.pairs.fasta_1.fasta $outfile_prs");
	}
	my $tmp_rename = "./" . "$prefix" . ".ALL.frags.ids.txt";
	runsys("mv ./ALL.frags.ids.txt $tmp_rename");
	$tmp_rename = "./" . "$prefix" . ".ALL.pairs.ids.txt";
	runsys("mv ./ALL.pairs.ids.txt $tmp_rename");
	$tmp_rename = "./" . "$prefix" . ".REDUCED_COVERAGE_ID_LIST.txt";
	runsys("mv ./REDUCED_COVERAGE_ID_LIST.txt $tmp_rename");
	
	
	if ((-s($outfile)) || ((-s($outfile_frg)) || (-s($outfile_prs))) ){
		if ($keep == 0){
			printl("Clean up...\n");
			runsys("rm INITBIN*");
			runsys("rm REDUCED_COVERAGE_ID_LIST.txt");
			runsys("rm uniq_counts_ids.txt");
			runsys("rm usable_reads_1.fasta");
			runsys("rm extractFasta.log");
			runsys("rm -r KMER_BINS");
			runsys("rm -r WHOLE_BINS");
			runsys("rm -r KMER_BINS_mp_frg");
			runsys("rm extractfasta.log");
		}else{
			runsys("mkdir INITBINS_CDHITEST");
			runsys("mv ./INITBIN* INITBINS_CDHITEST/");
		}
		if ($pairstatus==0){
			print "\n\nRun Successful! - See output:\n\t$outfile\n\nExiting.";
			print LOG "\n\nRun Successful! - See output:\n\t$outfile\n\nExiting.";
		}else{
			print "\n\nRun Successful! - See output:\n\t$outfile_frg\n\t$outfile_prs\n\nExiting.";
			print LOG "\n\nRun Successful! - See output:\n\t$outfile_frg\n\t$outfile_prs\n\nExiting.";
		}
	}else{
		printl("ERROR : NeatFreq encountered a problem and failed to produce an output file!  Investigate the end of log file NeatFreq.log.txt.\n\n");
	}
	
	#END OF MAIN:
	close LOG;
	
	my $newlog = "./" . "$prefix" . "_NeatFreq.LOG.txt";
	if (-s("$newlog")){
		system("rm $newlog");
	}
	system("mv ./NeatFreq.LOG.txt $newlog");
	
	print "Logs printed to $newlog\n\nEXITING.\n";
	exit;
}
