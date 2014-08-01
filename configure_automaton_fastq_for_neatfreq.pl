#!/usr/bin/perl -w
use strict;

sub outhelp(){
	print"FLAGS:\n-frg fragments.fastq\n-prs pairs.interleaved.fastq\n-out output_prefix\n-f QV offset\n-N = provide NeatFreq install location (if not updated in script)\n\nEXITING.\n";
	exit;
}

sub runsys {
	my $cmd = shift;
	print "++ RUN CMD : $cmd\n";
	system("$cmd");
}

MAIN : {
	print "CONFIGURE READS FOR NEATFREQ!\n";
	
	my $in; my $dir; my $out; my $dummy; my $prefix; my $infrg; my $inprs; my $offset=0; 
	my $NEW_NEATFREQ_INSTALL = `pwd`;
	chomp $NEW_NEATFREQ_INSTALL;
	while (@ARGV){
		$dummy=shift(@ARGV); 						#check -[a-z]
		if ( $dummy eq '-frg'){
			$infrg = shift(@ARGV);
			print "\nfrg in = $infrg\n";
		}elsif ( $dummy eq '-prs'){
			$inprs = shift(@ARGV);
			print "\nprs in = $inprs\n";
		}elsif ( $dummy eq '-out'){
			$prefix = shift(@ARGV);
			print "\nprefix in = $prefix\n";
		}elsif ( $dummy eq '-N'){
			$NEW_NEATFREQ_INSTALL = shift(@ARGV);
			print "NeatFreq install location in = $NEW_NEATFREQ_INSTALL\n";
			print "Forcing use of user-selected NeatFreq install directory : $NEW_NEATFREQ_INSTALL\n";
			if (!(-s("$NEW_NEATFREQ_INSTALL/NeatFreq.pl"))){
				print "\n\nERROR! : No NeatFreq install found in suggested install directory.  \n\nExiting...\n";
				exit;
			}
		}
		elsif ( $dummy eq '-f'){
			$offset = shift(@ARGV);
			print "\noffset in = $offset\n";
		}
		elsif ( $dummy eq '-h'){
			outhelp();
		}
		else{ 										#incorrect -[a-z] input
			print "*ERROR: Incorrect data entry $dummy.\n\n";
			outhelp();
		}
		if (!(-s("$NEW_NEATFREQ_INSTALL/NeatFreq.pl"))){
			print "\n\nERROR! : No NeatFreq install found in suggested install directory. ( $NEW_NEATFREQ_INSTALL/NeatFreq.pl )\n\nExiting...\n";
			exit;
		}
	}
	
	if ( (!($infrg)) || (!($inprs)) || (!($prefix)) ){
		print "COMMAND LINE INPUTS MISSING!\n\n";
		outhelp();
	}
	
	runsys("cat $infrg > all.Z.fastq");
	runsys("cat $inprs >> all.Z.fastq");
	if ($offset == 0){
		runsys("$NEW_NEATFREQ_INSTALL/lib/fastx_renamer -n COUNT -i all.Z.fastq -o all.COUNT.fastq");
	}else{
		runsys("$NEW_NEATFREQ_INSTALL/lib/fastx_renamer -n COUNT -Q $offset -i all.Z.fastq -o all.COUNT.fastq");
	}
	
	#my $runcmd = "grep -c '^" . "@" .  $infrg;
	#my $frag_count = `grep -c '^+' $infrg`;
	#my $frag_count = `/usr/local/packages/clc-ngs-cell/sequence_info $infrg | grep sequences | awk \'{print \$4}\'`;
	my $frag_count = `wc -l $infrg | awk \'{print \$1}\'`;
	chomp $frag_count;
	$frag_count = $frag_count/4;
	print "\n\nDEBUG : frag_count = $frag_count\n\n";
	
	#my $check_count = `grep -c '^+' all.COUNT.fastq`;
	#chomp $check_count;
	
	#if ((!(-s("all.COUNT.fastq"))) || ($check_count <= 1)){
	#	runsys("/usr/local/devel/BCIS/assembly/tools/multiLINEfasta_to_singleLINEfasta.pl all.Z.fasta > all.Z.1line.fasta");
	#	runsys("/usr/local/bin/fastx_renamer -n COUNT -i all.Z.1line.fasta -o all.COUNT.fasta");
	#}
		
		
	print "Extracting fragments ( $frag_count )...\n";
	
	my $frag_out = "$prefix" . ".frg.fastq";
	my $pair_out = "$prefix" . ".prs.fastq";
	open (GC, '< all.COUNT.fastq');
	open (FO, '> frag_out');
	open (PO, '> pairs_out');
	my $count=1;
	while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		if (($line =~ m/^@/m) && ($count <= ($frag_count*4))){
			print FO "$line\n";
		}elsif ($count <= ($frag_count*4)){
			print FO "$line\n";
		}else{
			print PO "$line\n";
		}
		$count++;
	}
	close PO;
	close FO;
	close GC;
	
	runsys("mv frag_out $frag_out");
	runsys("mv pairs_out $pair_out");
	runsys("rm all.Z.fastq");
	#runsys("rm all.Z.1line.fasta");
	runsys("rm all.COUNT.fastq");
}
