#!/usr/bin/perl -w

#[1] = input
#[2] = output

use strict;

################################################################################################################################################################################################################################################
# runsys : 
################################################################################################################################################################################################################################################
sub runsys {
	my $cmd = shift;
	print "++ RUN CMD : $cmd\n";
	system("$cmd");
}

MAIN : {
	`CURRENTPATH=\$PWD`;
	`COMMON_BIN=/usr/local/common`;
		
	#input variable declarations
	my ($reflen, $in, $out);
	
	
	#pull in vars from command line
	while (@ARGV){
		$in = shift(@ARGV);
		$out = shift(@ARGV);
	}

	if (-s("$out")){ system("rm $out"); }
	
	open (SIZE_IN, "< $in");
	open (SIZE_OUT, ">> $out");
	
	my $count = 0;
	my $prev_id;
	my $curr_id;
	
	while (defined(my $line = <SIZE_IN>)){ #get each line
		chomp $line;
		my @xs = split(/\s+/, $line);
		if ($count == 0){
			$curr_id = $xs[0];
			if ($xs[0] != 0){
				for (my $z = 0; $z < $xs[0]; $z++){
					print SIZE_OUT "$z\t+0\t1\tdeletethissequence!\n";
				}
				$prev_id = ($xs[0]-1);
			}else{
				$prev_id = 0;
			}
			print SIZE_OUT "$line\n";
		}else{
			$curr_id = $xs[0];
			if ($curr_id != $prev_id){ #end of sequence found
				if ($curr_id != ($prev_id+1)){ #if previous is not next sequentially...
					for (my $z = ($prev_id+1); $z < $curr_id; $z++){
						print SIZE_OUT "$z\t+0\t9999\tdeletethissequence!\n";
					}
				}
				print SIZE_OUT "$line\n";
			}else{
				print SIZE_OUT "$line\n";
			}
			$prev_id = $curr_id;
		}
		$count++;
		
	}	
	close SIZE_IN;
	close SIZE_OUT;
	

} # END MAIN
