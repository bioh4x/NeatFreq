#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;

my $inFasta = $ARGV[0];
my $baseName = basename($inFasta, qw/.fasta .fna/);
my $inQual = $baseName . ".qual";
my $outFastq = $baseName . ".fastq";

print "$inQual\t $inFasta \t $outFastq\n";
my %seqs;

$/ = ">";

open (FASTA, "<$inFasta");
my $junk = (<FASTA>);

my $fdef = "";
my $prev_fdef = "NULL";

while (my $frecord = <FASTA>) {
	chomp $frecord;
	
	my $prev_fdef = $fdef;
	my @seqLines = ();
	
	($fdef, @seqLines) = split /\n/, $frecord;
	my $seq = join '', @seqLines;
	
	if (($fdef eq $prev_fdef) && ($fdef ne "")){
		$fdef = "$fdef" . "_2";
	}
	
	$seqs{$fdef} = $seq;
}

close FASTA;

open (QUAL, "<$inQual");
$junk = <QUAL>;
open (FASTQ, ">$outFastq");

my $qdef = "";
my $prev_qdef = "NULL";

while (my $qrecord = <QUAL>) {
	chomp $qrecord;
	
	my $prev_qdef = $qdef;
	my @qualLines = ();
	
	($qdef, @qualLines) = split /\n/, $qrecord;
	my $qualString = join ' ', @qualLines;
	my @quals = split / /, $qualString;
	
	if (($qdef eq $prev_qdef) && ($qdef ne "")){
		$qdef = "$qdef" . "_2";
		print FASTQ "@","$qdef\n";
	}else{
		print FASTQ "@","$qdef\n";
	}
	
	print FASTQ "$seqs{$qdef}\n";
	print FASTQ "+\n";
	foreach my $qual (@quals) {
		print FASTQ chr($qual + 64);
	}
	print FASTQ "\n";
}

close QUAL;
close FASTQ;
