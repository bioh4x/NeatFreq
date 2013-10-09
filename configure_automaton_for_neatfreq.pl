#!/usr/bin/perl -w
use strict;

sub outhelp(){
	exit;
}

sub runsys {
	my $cmd = shift;
	print "++ RUN CMD : $cmd\n";
	system("$cmd");
}

MAIN : {
	print "CONFIGURE READS FOR NEATFREQ!\n";
	
	my $in; my $dir; my $out; my $dummy; my $prefix; my $infrg; my $inprs;
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
		}
		else{ 										#incorrect -[a-z] input
			print "*ERROR: Incorrect data entry $dummy.\n\n";
			outhelp();
		}
	}
	
	runsys("cat $infrg > all.Z.fasta");
	runsys("cat $inprs >> all.Z.fasta");
	runsys("/usr/local/bin/fastx_renamer -n COUNT -i all.Z.fasta -o all.COUNT.fasta");
	
	my $frag_count = `grep -c '>' $infrg`;
	chomp $frag_count;
	
	my $check_count = `grep -c '>' all.COUNT.fasta`;
	chomp $check_count;
	
	if ((!(-s("all.COUNT.fasta"))) || ($check_count <= 1)){
		runsys("/usr/local/devel/BCIS/assembly/tools/multiLINEfasta_to_singleLINEfasta.pl all.Z.fasta > all.Z.1line.fasta");
		runsys("/usr/local/bin/fastx_renamer -n COUNT -i all.Z.1line.fasta -o all.COUNT.fasta");
	}
		
	
		
	print "Extracting fragments ( $frag_count )...\n";
	
	my $frag_out = "$prefix" . ".frg.fasta";
	my $pair_out = "$prefix" . ".prs.fasta";
	open (GC, '< all.COUNT.fasta');
	open (FO, '> frag_out');
	open (PO, '> pairs_out');
	my $count=1;
	while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		if (($line =~ m/^>/m) && ($count <= ($frag_count*2))){
			print FO "$line\n";
		}elsif ($count <= ($frag_count*2)){
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
	runsys("rm all.Z.fasta");
	runsys("rm all.Z.1line.fasta");
	runsys("rm all.COUNT.fasta");
}