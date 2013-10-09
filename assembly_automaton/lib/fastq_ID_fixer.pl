#!/usr/bin/perl -w
use strict;

sub outhelp(){
	print "fastq_ID_fixer.pl -i input.fastq\nWhen you have an ID line but empty quality ID line... Fills in quality ID line.  Allows maq fastq extraction.\n";
	exit;
}

MAIN : {
	my $in; my $dir; my $out; my $dummy;
	while (@ARGV){
		$dummy=shift(@ARGV); 						#check -[a-z]
		if ( $dummy eq '-i'){
			$in = shift(@ARGV);
			#print "\nin = $in\n";
		}
		else{ 										#incorrect -[a-z] input
			print "*ERROR: Incorrect data entry $dummy.  Please see help:\n\n";
			outhelp();
		}
	}
	
	
	open (READ_IN, "$in");
	my $line;
	#print "Reading $in...\n";
	my $count = 1;
	my $id;
	while (defined(my $line = <READ_IN>)){ #get each line
		chomp $line;
		#if ($line[0] eq "@"){
		if ( $line =~ /^@/ ){
			$count++;
			print "$line\n";
			$id = substr($line, 1);
		}
		elsif( ($count == 3) && ( $line =~ /^\+/ ) ){
			my $tmpstr = "+" . "$id";
			print "$tmpstr\n";
			$count = 0;
		}
		else{
			$count++;
			print "$line\n";
		}
	}
}