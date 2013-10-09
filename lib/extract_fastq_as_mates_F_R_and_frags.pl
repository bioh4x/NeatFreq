#!/usr/bin/perl -w

# TAKES INPUT OF SORTED LIST OF IDS
# OUTPUTS INTERLEAVED PAIR LIST AND FRAGMENT LIST
# ANTICIPATES MATES WITH A /1 and /2 SUFFIX

use strict;

MAIN : {
`CURRENTPATH=\$PWD`;
`COMMON_BIN=/usr/local/common`;
	
#input variable declarations
my $names; my $prefix;

#pull in vars from command line
while (@ARGV){
	$names =shift(@ARGV); #check -[a-z]
	$prefix = shift(@ARGV);
}

my $cur_id = "null"; my $prev_id;
my $count = 0;
my $for_mates = "$prefix" . ".for_mate_ids.txt";
my $rev_mates = "$prefix" . ".rev_mate_ids.txt";
my $fragments = "$prefix" . ".fragment_ids.txt";
system("touch $for_mates");
system("touch $rev_mates");
system("touch $fragments");

open (SIZE_IN, "$names");
while (defined(my $line = <SIZE_IN>)){ #get each line
	chomp $line;
	my $X = $line;
	chop($X);
	chop($X);
	if ($cur_id ne "null"){
		#print "X: $X -- CUR_ID: $cur_id\n";
		if ($X eq "$cur_id"){
			$count = 0;
			system("echo $prev_id >> $for_mates");
			system("echo $line >> $rev_mates");
			
		}
		elsif( $count == 1){ #new id , mate expected
			system("echo $prev_id >> $fragments");
			$count = 1;
		}
		else{ #new id, no mate expected
			$count = 1;
		}
	}
	$prev_id = "$line";
	$cur_id = "$X";
}

} # END MAIN
