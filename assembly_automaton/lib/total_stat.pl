#!/usr/bin/perl -w

###################################################################
# total_stat.pl                                                   #
# 1/5/2011        Z. Goodwin                                      #
# version 2.0 - NOW UPDATED BY R. SANKA AND J. MCCORRISON!!!!     #
# rsanka@jcvi.org jmccorri@jcvi.org                               #
#                                                                 #
# Please see README  and the -h flag for more details             #
#                                                                 #
# Input       : contigs.fasta, scaffolds.fasta,                   #
#             : 454Reads.sff, mates.fastq/fasta,                  # 
#             : fragments.fastq/fasta                             #                 
#                                                                 #
# Output      : Tab-delimited file showing assembly stats         #
#                                                                 #
# Dependencies: CLC license.properties file in working dir        #
#             : clc_ref_assemble_long/short                       #
#             : split_sequences                                   #
#             : sizeFasta                                         #
#             : find_nX_genome_size                               #
#             : RealGC.pl                                         #
#             : SAMTools                                          #
#                                                                 #
###################################################################

# PACKAGES
#
use strict;
use Getopt::Long;
use File::Basename;
#use Config::IniFiles;  


##################### Definitions #####################

# CONSTANTS
my $DEFAULT_PERC_ID     = 0.95;
my $DEFAULT_PERC_LENGTH = 0.85;
my $DEFAULT_DIRECTION   = "ss";
my $DEFAULT_DIST_MODE   = "fb";
my $DEFAULT_COV_CUTOFF  = 5;

#time tracking contants;
my @prev_time;
my $time_bool=0;

## Software Locations ##
my $LICENSE_DIR           = "/usr/local/packages/clc-ngs-cell/license.properties";   
my $SAMTOOLS_DIR;
#update 8/18
if (-s("/usr/local/packages/samtools-0.1.18/bin/samtools")){
	print "CENTOS6!\n";
	$SAMTOOLS_DIR = "/usr/local/packages/samtools-0.1.18/bin/";
}elsif (-s("/usr/local/packages/samtools-0.1.13/bin/samtools")){
	print "CENTOS5!\n";
	$SAMTOOLS_DIR = "/usr/local/packages/samtools-0.1.13/bin/";
}else{
	printl("\n\nFAIL: No samtools versions found.\nExiting.\n\n"); 
	exit;
}
my $clc_ref_assemble_long = "/usr/local/packages/clc-ngs-cell/clc_ref_assemble_long";
my $castosam        = "/usr/local/packages/clc-ngs-cell/castosam";
my $assembly_info         = "/usr/local/packages/clc-ngs-cell/assembly_info";
my $assembly_table        = "/usr/local/packages/clc-ngs-cell/assembly_table";
my $sequence_info         = "/usr/local/packages/clc-ngs-cell/sequence_info";
my $split_sequences       = "/usr/local/packages/clc-ngs-cell/split_sequences";

# STANDARD OPTIONS
#
my ($prefix,$contigFile,$scaffFile,$expGenomeLength,$mapStatus,$help,$v,$noheader,$sfftitan,$sffflx,$pst_sample,$pst_project,$zero_x_ref,$insert_span_ref);
my $contigCount = 0;

# MATE-SPECIFIC OPTIONS 
#
my (%frg,%imates,%sffflx,%sfftitan,%frgFiles,%frmates,%ins,%mateOrient,%matePoints);
my (@frags,@fragsOnly,@mates,@reads,@ins,@dat) = ();

my $length = $DEFAULT_PERC_LENGTH;
my $identity = $DEFAULT_PERC_ID;
my $covCutoff = $DEFAULT_COV_CUTOFF;
my $lenCutoff = 0;

# GLOBALS
#
my ($keep,$autoLength) = (0) x 2; # flags for keeping the .cas file and using the contig span as genome length
my $OUTPUT = ""; 
my ($dummy,$ref,$opts,$version);
my $libID = 0;
my %READ_INFO;
my %MATEINFO;

# TIME STAMP
if (-e("./total_stat.time_stamp.LOG.txt")){ system("rm ./total_stat.time_stamp.LOG.txt"); }
system("touch ./total_stat.time_stamp.LOG.txt");
print "debug\n";
open LOG, ">> ./total_stat.time_stamp.LOG.txt";

# THE READFILE AND CLC-OPTIONS MENU
#
my $usage = "\ntotal_stat.pl -- VERSION 2  -- Generate basic statistics on a raw assembly".
"\n\nUSAGE:       perl total_stat.pl <options>\n\n".
"\nCOMMAND-LINE OPTIONS (must always be given):".
"\n\t-c <contigs.fasta>".
"\n\t-g [expected genome length | auto] Auto means that the cumulative contig span is used".
"\n\t             instead of the true expected genome length. Use with caution!\n".
"\n\nOTHER OPTIONS:".
"\n\t-p <filename prefix/organism name (no spaces)>".
"\n\t-s <scaffolds.fasta>".
"\n\t-h <print this help message>\n\n\n".
"READ MAPPING OPTIONS (provide only if mate pairs used for assembly):".
"\n\t-is         -is <lower_bound_mate_range> <upper_bound_mate_range>\n\t\t\t ALWAYS put this at the start of each read library".     
"\n\t-frg <frag1.fasta> <frag2.fasta> ...".
"\n\t-fr         -fr <forward.mates.fasta> <reverse.mates.fasta> (can also be in .fastq format)".
"\n\t-ii         <SINGLE, interleaved mate pair file in .fastq or .fasta format>".
"\n\t-sffflx     <454_FLX mates.sff>".
"\n\t-sfftitan   <454_TITANIUM mates.sff>".
#"\n\t-d          <mate pair orientation (can be ff,fb,bf or bb, default = ss)>".
#"\n\t-mp         <mate pair distance mode (can be ss,ee,es or se, default = fb)>".
"\n\t-cov <minimum coverage level> (for % of genome with >= X coverage stat, default = 5)".
"\n\t-l <percent read length to match (default = 0.85)>".
"\n\t-i <read percent identity (default = 0.95)>".
"\n\t-k (Keep .cas file from clc_ref_assemble_long reference assembly )".
"\n\t-v (print vertical output)".
"\n\t-nohead (do not show header in output)".
"\n\t-pst_project <project_id> (and print PST output, requires sample id)".
"\n\t-pst_sample <sample_id> (and print PST output, requires project id)".
"\n\t-zero_x <reference.fasta> (and print count of 0x bp)".
"\n\t-insert_span <reference.fasta> (and print average insert span, required input newbler-formatted mates in fasta format, warning: slow!)".
#NEW-JM
"\n\t-V    <Assembler Version (input by user, for db storage)>".
"\n\t-cutoff     <minimum contig length (default = 0, recommended = 500) (set to \"annot\" for annotation-ready contigs)>\n\n".

"For bug reports and updates, please see http://jira.jcvi.org/browse/BCIS-315\n";

#################### PARSE ALL INPUT FILES FROM CMD LINE ####################

print "### TOTAL_STAT.pl ###\n\n";

# Store command line options
my $retryCommand = join(" ",@ARGV);
my $logstr = "";

while (@ARGV){
    
    $dummy =shift(@ARGV); #check -[a-z]
    
    if(defined($dummy)){
		if ( $dummy eq '-c' ){ 
			
			$contigFile=shift(@ARGV); 
			print "\tUsing contig File : $contigFile\n";
			$logstr = $logstr."\tUsing contig File : $contigFile\n";
		}
		elsif ( $dummy eq '-p' ){
			
			$prefix=shift(@ARGV); 
			
			if($prefix =~ m/[!@#$%^&*()]/){   ### TEST
				print "\tInvalid file name prefix: $prefix";
			$logstr = $logstr."\tInvalid file name prefix: $prefix";
	    }else{
			print "\tUsing Prefix : $prefix\n";
			$logstr = $logstr."\tUsing Prefix : $prefix\n";
	    }
	}
	elsif ( $dummy eq '-g'){
	    
	    $expGenomeLength=shift(@ARGV); 
	    
	    # ADDED 7/7/2011: automatically run total_stat with 
	    # contig span as expected genome length (useful for large numbers of assemblies)
	    if( $expGenomeLength eq 'auto' ){
			$autoLength = 1;
			print "\n\t***\n\tWARNING: expected genome length has been set to 'auto.'".
		    "\n\tThe cumulative contig span may not represent the actual genome length!\n".
		    "\n\t***\n";
			$logstr = $logstr."\n\t***\n\tWARNING: expected genome length has been set to 'auto.'".
		    "\n\tThe cumulative contig span may not represent the actual genome length!\n".
		    "\n\t***";
	    }
		
	    print "\tExpected Genome Length : $expGenomeLength\n";
	    $logstr = $logstr."\tExpected Genome Length : $expGenomeLength\n";
	}
	elsif ( $dummy eq '-k'){$keep = 1;}   ## Set flag to keep reference assembly file (.cas file)
	elsif ( $dummy eq '-h'){die $usage;}  ## print help message
	elsif ( $dummy eq '-v'){$v=1;}
	elsif ( $dummy eq '-nohead'){$noheader=1}
	elsif ( $dummy eq '-l'){ 
	    
	    $length=shift(@ARGV);
	    print "Using % length : $length\n";
	    $logstr = $logstr."Using % length : $length\n";
	} 
	elsif ( $dummy eq '-i'){ 
		
	    $identity=shift(@ARGV);
	    print "Using % identity : $identity\n";
	    $logstr = $logstr."Using % identity : $identity\n";
	}
	elsif ( $dummy eq '-s'){ 
	    
	    $scaffFile=shift(@ARGV); 
	    print "\tUsing Scaffold File : $scaffFile\n";
	    $logstr = $logstr."\tUsing Scaffold File : $scaffFile\n";
	}
	elsif ( $dummy eq '-cov'  ){
	    $covCutoff=shift(@ARGV);
	    print "Using coverage cutoff : $covCutoff\n";
	    $logstr = $logstr."Using coverage cutoff : $covCutoff\n";
	}
	elsif ( $dummy eq '-cutoff'  ){
	    $lenCutoff=shift(@ARGV);
	    if ($lenCutoff eq "annot"){
			print "\nANNOTATION CHECK REQUESTED\t(500 bp, no terminal N's)\n\n";
			$logstr = $logstr."\nANNOTATION CHECK REQUESTED\t(500 bp, no terminal N's)\n\n";
	    }else{
	   		print "Using contig length curoff : $lenCutoff\n";
	    	$logstr = $logstr."Using contig length curoff : $lenCutoff\n";
		}
	}
	elsif ( $dummy eq '-V'  ){
	    $version=shift(@ARGV);
	    print "Using version number : $version\n";
	    $logstr = $logstr."Using version number : $version\n";
	}
	elsif ( $dummy eq '-is'){  # When this option appears, enable read mapping
	    
	    $mapStatus = 1;
	    
	    (@frags,@mates,@ins) = ();
	    
	    push(@ins,shift(@ARGV)); # Get the first insert length
	    push(@ins,shift(@ARGV)); # Get the second insert length
	    
	    print "\tFound Inserts.\n";
	    $logstr = $logstr."\tFound Inserts.\n";
		
	    chomp(@ins);
	    $libID++;
	    
	    # Initialize all the file name hashes for each library
	    $sfftitan{$libID}   ="";	
	    $sffflx{$libID}     ="";
	    $imates{$libID}     ="";
	    $ins{$libID}        ="";
	    $frgFiles{$libID}   ="";   
	    $frmates{$libID}    ="";
	    $mateOrient{$libID} ="";
	    $matePoints{$libID} ="";
	    
	    print "\n\t****\n\tInsert Length: @ins\n";
	    $logstr = $logstr."\n\t****\n\tInsert Length: @ins\n";
	}
	elsif ( $dummy eq '-frg'){      
	    $mapStatus = 1;
        do{
			$dummy=shift(@ARGV);       
			push(@fragsOnly, $dummy);
	    }until((!defined($dummy)) || ($dummy =~ m/^-/));
	    
	    unshift(@ARGV,$dummy);  # push the last option back onto ARGV so it can be read later
	    pop(@fragsOnly);        # pop the flag off the fragment array so it doesn't get read as an input file
		
	    if(defined($dummy)){
			my $pstr = "\tFound fragment files : ".join(",",@fragsOnly)."\n";		    
			print $pstr;		    
			$logstr = $logstr.$pstr;		    
	    }
	}
	elsif ( $dummy eq '-fr'){
	    
	    do{
			$dummy=shift(@ARGV);       
			push(@mates,$dummy);
	    }until((!defined($dummy)) || ($dummy =~ m/^-/));
		
	    unshift(@ARGV,$dummy);
	    pop(@mates);
		
	    if(defined($dummy)){
			my $pstr = "\tFound frmates : ".join(",",@mates)."\n";
			print $pstr;		
			$logstr = $logstr.$pstr;
	    }
	}
	elsif ( $dummy eq '-ii'){
	    $dummy=shift(@ARGV);       	    
	    $imates{$libID}=$dummy;
	    print "\tFound interleaved mates : $dummy\n";
	    $logstr = $logstr."\tFound interleaved mates : $dummy\n";
	    
	}
	elsif ( $dummy eq '-sffflx'){   
	    $dummy=shift(@ARGV);
	    $sffflx{$libID} = $dummy;
	    print "\tFound sff FLX : $dummy\n";
	    $logstr = $logstr."\tFound sff FLX : $dummy\n";
	}
	elsif ( $dummy eq '-sfftitan'){ 
	    $dummy=shift(@ARGV);
	    $sfftitan{$libID} = $dummy;
	    $logstr = $logstr."\tFound sff TITAN : $dummy\n";
	    print "\tFound sff TITAN : $dummy\n";
	}
	elsif ( $dummy eq '-d'){
	    $dummy=shift(@ARGV);        
	    $mateOrient{$libID}= $dummy;
	    $logstr = $logstr."\tFound Mate Orientation : $dummy\n";
	    print "\tFound Mate Orientation : $dummy\n";
	}
	elsif ( $dummy eq '-mp'){
	    $dummy=shift(@ARGV);       
	    $matePoints{$libID}=$dummy;
	    print "\tFound Mate Points : $dummy\n";
	    $logstr = $logstr."\tFound Mate Points : $dummy\n";
	}
	elsif( $dummy eq '-pst_project'){
		$pst_project=shift(@ARGV); 
		print "\tUsing PST Project ID : $pst_project\n";
		$logstr = $logstr."\tUsing PST Project ID : $pst_project\n";
	}
	elsif( $dummy eq '-pst_sample'){
		$pst_sample=shift(@ARGV); 
		print "\tUsing PST Sample ID : $pst_sample\n";
		$logstr = $logstr."\tUsing PST Sample ID : $pst_sample\n";
	}
	elsif( $dummy eq '-zero_x'){
		$zero_x_ref=shift(@ARGV); 
		print "\tUsing Zero X Comparison Reference : $zero_x_ref\n";
		$logstr = $logstr."\tUsing Zero X Comparison Reference : $zero_x_ref\n";
	}
	elsif( $dummy eq '-insert_span'){
		$insert_span_ref=shift(@ARGV); 
		print "\tUsing Insert Span Comparison Reference : $insert_span_ref\n";
		$logstr = $logstr."\tUsing Insert Span Comparison Reference : $insert_span_ref\n";
	}
	#-zero_x <reference.fasta> (and print count of 0x bp)".
	#"\n\t-insert_span
	else{
	    $logstr = $logstr."Unrecognized option : $dummy\n";   ### TEST
	    die "Unrecognized option : $dummy\n";
	}
}

if($libID != 0){
	$ins{$libID}      =join(",",@ins);
	$frmates{$libID}  =join(",",@mates);
}
}#End while




# NEW - JM - 8/24/11
# CHECK FOR CUTOFF, RUN NOW!
if ($lenCutoff ne 0){
	if ($lenCutoff eq "annot"){
		print "\nANNOTATION CHECK : Removing contigs at less than 500 length and deleting terminal N's.\n";
		
		#TIMESTAMP
		supply_time();
		
		$contigFile = runANNOTcheck($contigFile);
	}else{
		print "\nRemoving contigs at less than $lenCutoff length.\n";
		
		#TIMESTAMP
		supply_time();
		
		$contigFile = runLENCUTOFF($contigFile, $lenCutoff);
		
	}
}


# Wrap up the fragment files that don't belong to an insert library
#
$frgFiles{'FRAGSONLY'} = join(",",@fragsOnly);

print "\tDONE - input files found\n\t****\n";
$logstr = $logstr."\tDONE - input files found\n\t****\n";

%READ_INFO = ('insert'      => \%ins,
'frags'       => \%frgFiles,
'frmates'     => \%frmates,
'sffflx'      => \%sffflx,
'sfftitan'    => \%sfftitan,
'interleaved' => \%imates,
'orient'      => \%mateOrient,
'distmode'    => \%matePoints
);


############################ MAIN ############################

my $tmp_uptime = localtime();

#TIMESTAMP
supply_time();

# Initialize all stats.  ("-" for null value ... implies that the stat was not run)
my %STATS = ('Run Prefix' => "-",
'Contig File' => "-",
'Genome Length' => "-",
'Contig count'  => "-",
'Contig n50' => "-",
'Contig n75' => "-",
'Contig n90' => "-",
'Len. of longest contig' => "-",
'Total Contig Len.' => "-",
'Average Contig Len.' => "-",
'Genome %GC Content' => "-",
'Scaffold File' => "-",
'Scaffold count' => "-",
'Scaffold n50' => "-",
'Scaffold n75' => "-",
'Scaffold n90' => "-",
'Intra-scaffold gaps' => "-",
'Gaps per 5kb' => "-",
'Num. of bases in reads' => "-",
'Num. of reads used' => "-",
'Avg. read length' => "-",
'Total Mates' => "-",
'Average genome coverage' => "-",
'Percent coverage' => "-",
'Missing Reference BPs' => "-",
'Average Insert Span' => "-");

#'Percent coverage' => "-");
# NEW
# PROCESS CONTIG FILE

#
if(($contigFile) && (-s $contigFile)){
    
    # Use contig file name as default prefix if prefix is not specified.
    # 
    if(!$prefix){$prefix = basename($contigFile);}
    
    # Print all found options to log, once the prefix is known.
    #
    printlog($logstr);
	
    # PREPARE OUTPUT AND LOGFILE
    #
	#    if((!$v) && (!$noheader)){ print_header($contigFile,$scaffFile,$mapStatus,$prefix,$covCutoff); }
	#    if(-e "./$prefix.log"){runsys("/bin/rm ./$prefix.log");}    
    
    printlog("Using default filename prefix: $contigFile\n");
    printlog("Using Contig File: ".$contigFile."\n");
    
    # If expected genome length is not given, calculate predicted genome length and end the program
    # Users may use the predicted genome length at their own risk.  
    #
    if(!$expGenomeLength){
		
		# Get expected genome length
		my $tmp = seqInfo($contigFile);
		@dat = @{$tmp};
		$expGenomeLength = $dat[0];
		chomp($expGenomeLength);
		
		# Create a new command with the expected genome size for easy re-run
		my @CMD = split(/\s+/,$retryCommand);
		unshift(@CMD,$expGenomeLength);
		unshift(@CMD,"-g");
		$retryCommand = join(" ",@CMD);
		
		die "\nExpected genome length not specified.\n".
	    "Predicted genome length = $expGenomeLength\n".
	    "Try re-running total_stat with the following command:".
	    "total_stat.pl $retryCommand\n";
		
    }
	
    # if '-g auto' is specified, use the cumulative contig
    # span as the expected genome length
    if($autoLength){
		
		my $tmp = seqInfo($contigFile);
		@dat = @{$tmp};
		$expGenomeLength = $dat[0];
		chomp($expGenomeLength);
    }
    
    # Get Contig Stats
    #
    $contigCount = calculateContigStats($contigFile,$expGenomeLength,$v,$prefix);
    
}else{ die "\nERROR : Contig file not found.\n"; }

# PROCESS SCAFFOLD FILE
printl("Scaffold stats\n");
supply_time();
if($scaffFile){ 
    
    if((-s $scaffFile) && (-e $scaffFile)){
		
		printlog("Using Scaffold File: ".$scaffFile."\n");
		my $scaffoldCount = calculateScaffoldStats($scaffFile,$contigCount,$expGenomeLength,$v);
		
    }else{
		
		printlog("\nScaffold file empty or non-exsitent.  Please check command and file scaffold file contents.\n\n");
		die "\nScaffold file empty or non-exsitent.  Please check command and file scaffold file contents.\n\n";
		
    }
}

# PROCESS MATE FILES
printl("PROCESS MATE FILES\n");

#TIMESTAMP
supply_time();

if( ($libID != 0)  || ($READ_INFO{'frags'}{'FRAGSONLY'} ne "") ){
    
    # Check if percent identity, percent length, mate direction and distance mode are specified
    # If not, use default values.
    #
    if(!$identity){
		printlog("Percent identity not specified.  Using default.\n");
		$identity = $DEFAULT_PERC_ID;
    }
    if(!$length){
		printlog("Percent length not specified.  Using default.\n");
		$length = $DEFAULT_PERC_LENGTH;
    }
    
    # Check for presence of license.properties, can't run the clc mapper without it!
    #
    unless (-e "./license.properties"){
		
		printlog("License file missing from working dir.  Copying license file to current directory...  \n");
		
		# copy the license file to the working directory
		if (-e $LICENSE_DIR){ `cp $LICENSE_DIR ./`; }
		else{
			printlog("\nLicense file not present in ".$LICENSE_DIR."\n\n");
			die "\nLicense file not present in ".$LICENSE_DIR."\n\n";
		}
    }
    
	## IMPORTANT!! ##
	#    Convert 454 files, get read statistics, map reads, then calculate coverage stats
	#   
	my $ref=checkFiles(\%READ_INFO);    
    
	$ref=parse454(\%READ_INFO,$prefix);
	%READ_INFO = %{$ref};    
	
	printl("GetReadStats\n");
	supply_time();
    getReadStats(\%READ_INFO,$v);
    printl("runCLCMapping\n");
    supply_time();
    runCLCMapping($contigFile,$length,$identity,\%READ_INFO,$prefix);
    printl("coverageStats\n");
    supply_time();
    coverageStats($contigFile,$prefix,$v,$covCutoff);
    printl("coverageStats complete\n");
    supply_time();
    
    if ($zero_x_ref){
    	printl("ZeroXAnalysis\n");
    	supply_time();
    	my $tmpzzz = ZeroXAnalysis($contigFile,$prefix,$zero_x_ref);
    	printl("coverageStats complete\n");
    	$STATS{'Missing Reference BPs'} = "$tmpzzz";
    	supply_time();
    }
#	if ($insert_span_ref){
#		printl("InsertSpanCheck\n");
#    	supply_time();
#    	InsertSpanCheck($contigFile,$prefix,$v,$covCutoff);
#    	printl("coverageStats complete\n");
#	}

}

# PRINT ALL STATISTICS TO SCREEN
printStats(\%STATS,$contigFile,$scaffFile,$mapStatus,$covCutoff,$prefix,$v,$noheader,$pst_sample,$pst_project);
my $tmp_end_time = localtime();
printlog("\nSTART TIME: ".$tmp_uptime."\n"."END  TIME: ".$tmp_end_time."\n");

# Cleanup
#if(-e "./$contigFile.fai"){system "/bin/rm ./$contigFile.fai";}
#if(-e "./$prefix.fai"){system "/bin/rm ./$prefix.fai";}

# put samtools error into log
if(-e "./samtools.err"){
    printlog("\n### SAMTOOLS ERRORS ###\n");
    my $out = `/bin/cat ./samtools.err`;
    printlog($out);
    printlog("\n#######################\n");
	#    runsys("/bin/rm ./samtools.err");
	

}


######################### END OF MAIN ############################

######################### SUBROUTINES ############################



#cutoff of contig length
sub runANNOTcheck{
	$contigFile = shift;
	
	runsys("/usr/local/devel/BCIS/assembly/tools/clean_contig_ends_fasta.pl -i $contigFile -o tmp.clean.fasta");
	
	runsys("/usr/local/devel/BCIS/assembly/tools/extractFasta -i tmp.clean.fasta -minsize 500 -o annotation_ready.ctg.fasta");
	
	runsys("mv annotation_ready.ctg.fasta_1.fasta annotation_ready.ctg.fasta");
	
	printl("++ RUN CMD : grep -c '>' annotation_ready.ctg.fasta\n");
	my $count = `grep -c '>' annotation_ready.ctg.fasta\n`;
	chomp $count;
	print "$count records larger than 500 bp written to file.\n";
	
	#TIMESTAMP
	printl("END OF runANNOTcheck\n");
	supply_time();
	
	return "annotation_ready.ctg.fasta";
}


#cutoff of contig length
sub runLENCUTOFF{
	$contigFile = shift;
	$lenCutoff = shift;
	
	my $SIZEFASTA_LOCATION = "/usr/local/devel/BCIS/assembly/tools";
	my $EXTRACTFASTA_LOCATION = "/usr/local/devel/SE/bin";
	
	my $count = 0;
	
	#runsys("cat $contigFile > $out_dir/$pre2.ctgANDdeg.fasta");
	printl("++ RUN CMD : grep -c '>' $contigFile\n");
	$count = `grep -c '>' $contigFile`;
	chomp $count;
	print "$count contig records found in original file.\n";
	
	#Parse out only those which are over 500 bp.
	runsys("$SIZEFASTA_LOCATION/sizeFasta $contigFile > all.lst");
	#create new output file
	open (FPLUS, "> 500plus.lst");
	open (SIZE_IN, "all.lst");
	while (defined(my $line = <SIZE_IN>)){ #get each line
		chomp $line;
		my @xs = split(/\s+/, $line);
		if ($xs[1] >= 500){ print FPLUS "$xs[0]\n";}
	}
	#clean up
	#		runsys("rm all.lst");
	close FPLUS;
	
	#split off of fasta file
	runsys("$EXTRACTFASTA_LOCATION/extractFasta -i $contigFile -idlist 500plus.lst -o z");
	runsys("mv z_1.fasta out.500bpContigs.fasta");
	$count = `grep -c '>' out.500bpContigs.fasta\n`;
	chomp $count;
	print "$count records larger than 500 bp written to file.\n";
	
	#TIMESTAMP
	printl("END OF runLENCUTOFF\n");
	supply_time();
	
	return "out.500bpContigs.fasta";

}



# Purpose : Check the if each input file in each library exists.  If not, kill the program
# Input   : %READ_INFO (the hash containing all read file paths)
sub checkFiles{
    my $ref   = shift;
    my %READS = %{$ref};
	
    my %frags    = %{$READS{'frags'}};
    my %frmates  = %{$READS{'frmates'}};  # Separate forward and reverse mate files
    my %sffflx   = %{$READS{'sffflx'}};
    my %sfftitan = %{$READS{'sfftitan'}};
    my %ileaved  = %{$READS{'interleaved'}};
    my @fragsOnly = split(/\,/,$frags{'FRAGSONLY'});
    
    my @err = (); # Missing files
    my @emp = (); # Empty files
	
    printlog("Now checking read files...\n");
	
    foreach my $fragsOnlyFile(@fragsOnly){
        if($fragsOnlyFile){
            if(-e $fragsOnlyFile){ printlog("Found fragment-only file $fragsOnlyFile \n");}
            else{push(@err,$fragsOnlyFile);}
            
            if(-z $fragsOnlyFile){push(@emp,$fragsOnlyFile);}
        }
    }
	
    for(my $i=1;$i<=$libID;$i++){
		
		my (@frags,@frmates,@ileaved,@tmp,@sffflx,@sfftitan) = ();
		
		@frags     = split(/\,/,$frags{$i});
		@frmates   = split(/\,/,$frmates{$i});
		@ileaved   = split(/\,/,$ileaved{$i});
		@sffflx    = split(/\,/,$sffflx{$i});
		@sfftitan  = split(/\,/,$sfftitan{$i});
		
		foreach my $frg(@frags){
			if($frg){
				if(-e $frg){ printlog("Found fragment-from-paired file $frg \n");}
				else{push(@err,$frg);}
				
				if(-z $frg){push(@emp,$frg);}
			}
		}
		
		@tmp = (@frmates,@ileaved,@sffflx,@sfftitan);
		
		foreach my $f(@tmp){
			if($f){
				if(-e $f){ printlog("Found paired file $f \n");}
				else{push(@err,$f);}
				
				if(-z $f){push(@emp,$f);}
			}
		}
    }
	
    # Show the missing files
    if(scalar(@err) > 0){
		foreach my $e(@err){
			print STDERR "File not found : $e\n";
			printlog("File not found : $e\n");
		}
		die "\nExiting.\n";
    }
	
    # Show the empty files
    if(scalar(@emp) > 0){
		foreach my $m(@emp){
			printlog("File not found : $m\n");
			print STDERR "Empty file : $m\n";
		}
		die "\nExiting.\n";
    }
	
	
	
    return;
}

# Purpose : Calculate the nX value of a set of sequences.  
# Input   : sequence file, genome size and the % cutoff for the nX
#
sub getNX{
	
    my $file = shift;
    my $x    = shift;
    my $genomeSize = shift;
    my @trackCTG;
	
    # Extract relevant fields from from sequence_info and store them in @out
    	printl("++ RUN CMD : $sequence_info -l $file\n");
    my @out = `$sequence_info -l $file`;    
    @out = @out[13..$#out];
	
    foreach my $l(@out){
		chomp($l);
		my @dat = split(/\s+/,$l);
		my $tmp = $dat[$#dat];
		if(defined($tmp)){push(@trackCTG,$tmp);}
    }
    
    my $find_n50=0;
    foreach my $size(sort {$b<=>$a} @trackCTG) {
		$find_n50 += $size;
		
		if($find_n50 >= (($genomeSize/100)* $x) ) {
			return $size;
		}
    }
    
    # Fail condition -- return this if expected genome size is larger than 50%/75%/90% of assembly length
    return "-";
}

# Purpose : This subroutine prints a line to an output file
# Input   : String to print to command line
#
sub printlog{
	
    my $item = shift;
    open(LOG2,">>./$prefix.log") or die "Couldn't open the log file!\n";
    print LOG2 "$item";
    close(LOG2);
	
}

# Process SFF Formatted files
#
sub parse454{
    my $ref  = shift;
    my $pre  = shift;
    my %info = %{$ref};
	
    my %fragments     = %{$info{'frags'}};
    my %frmates       = %{$info{'frmates'}};
    my %interleaved   = %{$info{'interleaved'}};
    my %sffflx        = %{$info{'sffflx'}};
    my %sfftitan      = %{$info{'sfftitan'}};
	
    for(my $i=1;$i<=$libID;$i++){
        #my @frags   = split(/\,/,$fragments{$i});
        my $titan   = $sfftitan{$i};
        my $flx     = $sffflx{$i};
        #my @ileaved = split(/\,/,$interleaved{$i});
		
	    my $linker;
        my $sffFile;
        my $fragname = "$pre.$i.frags.fasta";
        my $matename = "$pre.$i.mates.fasta";
        
        if($flx) {
            if((-e $flx) && (-s $flx)) {
                $linker = "flx";
                $sffFile = $flx;
	        }
            else {
                print "\nWARNING:  $flx was stated on the command line but file is not present. Please check command and file contents.\n\n";
                printlog("\nWARNING:  $flx was stated on the command line but file is not present. Please check command and file contents.\n\n");
	        }
	    }
	    elsif($titan){
	        if((-e $titan) && (-s $titan)){
		        $linker = "ti";
                $sffFile = $titan;
	        }
            else {
		        print "\nWARNING:  $titan was stated on the command line but file is not present. Please check command and file contents.\n\n";
		        printlog("\nWARNING:  $titan was stated on the command line but file is not present. Please check command and file contents.\n\n");
	        }
	    }
        
        # If there is a linker and sffFile ready, do the following:
        if (($linker) && ($sffFile)) {
            
            # Convert sff files to reads and fragments
                	printl("++ RUN CMD : $split_sequences -i $sffFile -s $fragname -p $matename -d $linker -m 25\n");
            `$split_sequences -i $sffFile -s $fragname -p $matename -d $linker -m 25`;
		    printlog("SPLIT_SEQUENCES COMMAND: $split_sequences -i $sffFile -s $fragname -p $matename -d $linker -m 25\n");
            
            # Check if mate and fragment parsing was successful, if not, user should check the files
        	if(( !(-e $fragname)) && ( !(-e $matename))){
		        printlog("\tFailed to parse sff file for mates and fragments.\n");
		        die "\tFailed to parse sff file for mates and fragments.\n";
		    }
			
    	    if((-e $fragname) && (-s $fragname)){
	            $fragments{$i} = $fragname;
    	    }
			
    	    if((-e $matename) && (-s $matename)){
	            $interleaved{$i} = $matename;
    	    }
        }
    }
    
    # Rebuild the read files hash
    $info{'interleaved'} = \%interleaved;
    $info{'frags'} = \%fragments;
    
    return \%info;
}


sub printStats{
	
    my $ref = shift;
    my %STATISTICS = %{$ref};
    my $ctg = shift;
    my $scf = shift;
    my $map = shift;
    my $covCutoff = shift;
    my $pre = shift;
    my $v   = shift;
    my $nohead = shift;
    my $pst_sample = shift;
    my $pst_project = shift;
    my $statStr = "";
    
    if (($pst_sample) && ($pst_project)){
    	printl("printStats - PST :\n");
    
    	$statStr = "$pst_sample, $pst_project, Run Prefix, ".$STATISTICS{"Run Prefix"}.
			"\n$pst_sample, $pst_project, Contig File, ".$STATISTICS{"Contig File"}.
			"\n$pst_sample, $pst_project, Genome Length, ".$STATISTICS{"Genome Length"}.
			"\n$pst_sample, $pst_project, Contig count, ".$STATISTICS{"Contig count"}.
			"\n$pst_sample, $pst_project, Contig n50, ".$STATISTICS{"Contig n50"}.
			"\n$pst_sample, $pst_project, Contig n75, ".$STATISTICS{"Contig n75"}.
			"\n$pst_sample, $pst_project, Contig n90, ".$STATISTICS{"Contig n90"}.
			"\n$pst_sample, $pst_project, Genome %GC Content, ".$STATISTICS{"Genome %GC Content"}.
			"\n$pst_sample, $pst_project, Len. of longest contig, ".$STATISTICS{"Len. of longest contig"}.
			"\n$pst_sample, $pst_project, Total Contig Len., ".$STATISTICS{"Total Contig Len."}.
			"\n$pst_sample, $pst_project, Average Contig Len., ".$STATISTICS{"Average Contig Len."};
			
			
			$statStr = $statStr."\n$pst_sample, $pst_project, Scaffold File, ".$STATISTICS{"Scaffold File"}.
			"\n$pst_sample, $pst_project, Scaffold count, ".$STATISTICS{"Scaffold count"}.
			"\n$pst_sample, $pst_project, Scaffold n50, ".$STATISTICS{"Scaffold n50"}.
			"\n$pst_sample, $pst_project, Scaffold n75, ".$STATISTICS{"Scaffold n75"}.
			"\n$pst_sample, $pst_project, Scaffold n90, ".$STATISTICS{"Scaffold n90"}.
			"\n$pst_sample, $pst_project, Intra-scaffold gaps, ".$STATISTICS{"Intra-scaffold gaps"}.
			"\n$pst_sample, $pst_project, Gaps per 5kb, ".$STATISTICS{"Gaps per 5kb"};
			
			$statStr = $statStr."\n$pst_sample, $pst_project, Num. of bases in reads, ".$STATISTICS{"Num. of bases in reads"}.
			"\n$pst_sample, $pst_project, Num. of reads used, ".$STATISTICS{"Num. of reads used"}.
			"\n$pst_sample, $pst_project, Avg. read length, ".$STATISTICS{"Avg. read length"}.
			"\n$pst_sample, $pst_project, Total Mates, ".$STATISTICS{"Total Mates"}.
			"\n$pst_sample, $pst_project, Average genome coverage, ".$STATISTICS{"Average genome coverage"}.
			"\n$pst_sample, $pst_project, Percent of genome with >= $covCutoff coverage, ".$STATISTICS{"Percent coverage"}."\n";
			#NEW-JM
			if ($version){
				$statStr = $statStr."$pst_sample, $pst_project, Assembler Version, " . "$version";
			}
			if ($STATISTICS{"Missing Reference BPs"} ne "-"){
				my $tmpvv = $STATISTICS{"Missing Reference BPs"};
				$statStr = $statStr."$pst_sample, $pst_project, Perc. of Reference BPs Missing, $tmpvv";
			}
			
		open(PST_OUT,">./$pre.PST.csv") or die "Could not open ./$pre.PST.tsv!";    
    	print $statStr;
    	print "\n";
    	print PST_OUT $statStr;
    	print PST_OUT "\n";
    	close(PST_OUT);
    }
    elsif( ($pst_sample) || ($pst_project) ){
    	printl("Error on PST print : only project or sample supplied, but not both!\nExiting...\n");
    	exit;
    }
    	
    if($v){
        #newline-seperated
		if($ctg){
			printl("printStats - VERTICAL :\n");
			$statStr = "Run Prefix\t".$STATISTICS{"Run Prefix"}.
			"\nContig File\t".$STATISTICS{"Contig File"}.
			"\nGenome Length\t".$STATISTICS{"Genome Length"}.
			"\nContig count\t".$STATISTICS{"Contig count"}.
			"\nContig n50\t".$STATISTICS{"Contig n50"}.
			"\nContig n75\t".$STATISTICS{"Contig n75"}.
			"\nContig n90\t".$STATISTICS{"Contig n90"}.
			"\nGenome %GC Content\t".$STATISTICS{"Genome %GC Content"}.
			"\nLen. of longest contig\t".$STATISTICS{"Len. of longest contig"}.
			"\nTotal Contig Len.\t".$STATISTICS{"Total Contig Len."}.
			"\nAverage Contig Len.\t".$STATISTICS{"Average Contig Len."};
			
			
			$statStr = $statStr."\nScaffold File\t".$STATISTICS{"Scaffold File"}.
			"\nScaffold count\t".$STATISTICS{"Scaffold count"}.
			"\nScaffold n50\t".$STATISTICS{"Scaffold n50"}.
			"\nScaffold n75\t".$STATISTICS{"Scaffold n75"}.
			"\nScaffold n90\t".$STATISTICS{"Scaffold n90"}.
			"\nIntra-scaffold gaps\t".$STATISTICS{"Intra-scaffold gaps"}.
			"\nGaps per 5kb\t".$STATISTICS{"Gaps per 5kb"};
			
			$statStr = $statStr."\nNum. of bases in reads\t".$STATISTICS{"Num. of bases in reads"}.
			"\nNum. of reads used\t".$STATISTICS{"Num. of reads used"}.
			"\nAvg. read length\t".$STATISTICS{"Avg. read length"}.
			"\nTotal Mates\t".$STATISTICS{"Total Mates"}.
			"\nAverage genome coverage\t".$STATISTICS{"Average genome coverage"}.
			"\nPercent of genome with >= $covCutoff coverage\t".$STATISTICS{"Percent coverage"}."\n";
			#NEW-JM
			if ($version){
				$statStr = $statStr."Assembler Version\t" . "$version";
			}
			if ($STATISTICS{"Perc. of Missing Reference BPs"} ne "-"){
				my $tmpvv = $STATISTICS{"Missing Reference BPs"};
				$statStr = $statStr."Perc. of Reference BPs Missing\t$tmpvv";
			}
				
		}
    }elsif(!$v){
    	printl("printStats - HORIZONTAL :\n");
		if(!$nohead){
			if($ctg){
				$statStr = "Run Prefix\t".
				"Contig File\t".
				"Genome Length\t".
				"Contig count\t".
				"Contig n50\t".
				"Contig n75\t".
				"Contig n90\t".
				"Len. of longest contig\t".
				"Total Contig Len.\t".
				"Average Contig Len.\t".
				"Genome %GC Content\t";
				
				$statStr = $statStr."Scaffold File\t".
				"Scaffold count\t".
				"Scaffold n50\t".
				"Scaffold n75\t".
				"Scaffold n90\t".
				"Intra-scaffold gaps\t".
				"Gaps per 5kb\t";
				
				$statStr = $statStr."Num. of bases in reads\t".
				"Num. of reads used\t".
				"Avg. read length\t".
				"Total Mates\t".
				"Average genome coverage\t".
				"Percent of genome with >= $covCutoff X coverage";
				#NEW-JM
				if ($version){
					$statStr = $statStr."\tAssembler Version";
				}
				if ($STATISTICS{"Missing Reference BPs"} ne "-"){
					$statStr = $statStr."\tPerc. of Reference BPs Missing";
				}
			}
		}
		
		if($ctg){
			$statStr = $statStr."\n".$STATISTICS{"Run Prefix"}."\t".
			$STATISTICS{"Contig File"}."\t".
			$STATISTICS{"Genome Length"}."\t".
			$STATISTICS{"Contig count"}."\t".
			$STATISTICS{"Contig n50"}."\t".
			$STATISTICS{"Contig n75"}."\t".
			$STATISTICS{"Contig n90"}."\t".
			$STATISTICS{"Len. of longest contig"}."\t".
			$STATISTICS{"Total Contig Len."}."\t".
			$STATISTICS{"Average Contig Len."}."\t".
			$STATISTICS{"Genome %GC Content"}."\t";
			
			$statStr = $statStr.$STATISTICS{"Scaffold File"}."\t".
			$STATISTICS{"Scaffold count"}."\t".
			$STATISTICS{"Scaffold n50"}."\t".
			$STATISTICS{"Scaffold n75"}."\t".
			$STATISTICS{"Scaffold n90"}."\t".
			$STATISTICS{"Intra-scaffold gaps"}."\t".
			$STATISTICS{"Gaps per 5kb"}."\t";
			
			$statStr = $statStr.$STATISTICS{"Num. of bases in reads"}."\t".
			$STATISTICS{"Num. of reads used"}."\t".
			$STATISTICS{"Avg. read length"}."\t".
			$STATISTICS{"Total Mates"}."\t".
			$STATISTICS{"Average genome coverage"}."\t".
			$STATISTICS{"Percent coverage"};
			#NEW-JM
			if ($version){
				$statStr = $statStr."\t$version";
			}
			if ($STATISTICS{"Missing Reference BPs"} ne "-"){
				my $tmpvv = $STATISTICS{"Missing Reference BPs"};
				$statStr = $statStr."\t$tmpvv";
			}
		}
    }
    
    
    # PRINT TO .TSV FILE AND SCREEN
    open(OUT,">./$pre.tsv") or die "Could not open ./$pre.tsv!";    
    print $statStr;
    print "\n";
    print OUT $statStr;
    print OUT "\n";
    close(OUT);
}


# Purpose : Print a header to output based on options given
# Input   : Contig file name, scaffold file name, mapping status
# 
sub print_header{
 
	 my $contigFile = shift;
	 my $scaffFile  = shift;
	 my $mapStatus  = shift;
	 my $prefix     = shift;
	 
	 if(($contigFile) && (-s $contigFile)){
	 $OUTPUT = "Run Name".
	 "\tContig File".
	 "\tGenome Length".
	 "\tContig count".
	 "\tContig_n50".
	 "\tContig_N75".
	 "\tContig_n90".
	 "\tGC Content".
	 "\tLength of longest contig".
	 "\tTotal contig len.".
	 "\tAvg. contig len.";
	 }
	 
	 if(($scaffFile) && (-s $scaffFile)){
	 $OUTPUT = "$OUTPUT".
	 "\tScaffold File".
	 "\tScaffold count".
	 "\tScaffold_n50".
	 "\tScaffold_n75".
	 "\tScaffold_n90".
	 "\tIntra-scaffold Gaps".
	 "\tGaps per 5kb";
	 }
	 
	 if($mapStatus){
	 $OUTPUT = "$OUTPUT".
	 "\tTotal num. of bases input".
	 "\tTotal Num. of reads input".
	 "\tAvg. read len\tTotal Mates".
	 "\tAvg. contig coverage\t".
	 "\% of genome with >= $covCutoff reads supporting each base";
	 }
	 
	 $OUTPUT = "$OUTPUT"."\n";
 }

# Purpose : Calculate contig n50, n75, n90, largest contig, GC content, total and average contig length
# Input   : $contigFile, $genomeLength (user-defined), vertical input flag 
# 
sub calculateContigStats{
	
    my $contigFile    = shift;
    my $genomeLength  = shift;
    my $v             = shift;
    my $pre           = shift;
    my ($contigCount,$avgCtgLen,$totalCtgLen,$largestContig,$GCcontent,$GC,$contig_n50,$contig_n75,$contig_n90,$G,$C) = (0) x 12;
	
    # Get all sequence info using clc's sequence_info script
    my $tmp = seqInfo($contigFile);
    my @dat = @{$tmp};
	
    # Stats come out of seqInfo in this order: (totalLength,maxCtgLength,minCtgLength,avgCtgLength)
    $totalCtgLen   =  $dat[0];
    $largestContig =  $dat[1];
    $avgCtgLen     =  $dat[3];
    $contigCount   =  $dat[4]; 
	
	# Check if contig span is less than 50% of expected genome length
	#
    if($totalCtgLen < ($genomeLength/2)){
		print "WARNING:  Total contig span is less than 50% of the total genome length. Using contig span as expected genome length\n";
		printlog("WARNING:  Total contig span is less than 50% of the total genome length.  Using contig span as expected genome length\n");
		$genomeLength = $totalCtgLen;
    }
	
	
	# Get nX Values
	#
    $contig_n50    = getNX($contigFile,50,$genomeLength);
    $contig_n75    = getNX($contigFile,75,$genomeLength);
    $contig_n90    = getNX($contigFile,90,$genomeLength);
	
	# Calculate GC content
	#
    $G = 0;
    $C = 0;
    open(GC,"<$contigFile");
    while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		if ($line =~m/^>/){
			#               print "XXXXXXXXXXXXXXXXXXXXXXXX\n";
		}
		else{
			
			$C += ($line =~ tr/C//);
			$G += ($line =~ tr/G//);
		}
    }
    close(GC);
	
    $GC = $G + $C;
    $GCcontent = ($GC/$totalCtgLen)*100;
    $GCcontent = sprintf("%.2f", $GCcontent);
    
	# Format output from find_nX
	#
    chomp($contig_n50);
    chomp($contig_n75);
    chomp($contig_n90);
    chomp($contigCount);
    chomp($largestContig);
    chomp($GCcontent);
    chomp($totalCtgLen);
	
	# Format Output -- if vertical option specified, print vertically.
	# if not, print horizontally
	#
    $STATS{"Run Prefix"} = $pre;
    $STATS{"Contig File"} = $contigFile;
    $STATS{"Genome Length"} = $genomeLength;
    $STATS{"Contig count"} = $contigCount;
    $STATS{"Contig n50"} = $contig_n50;
    $STATS{"Contig n75"} = $contig_n75;
    $STATS{"Contig n90"} = $contig_n90;
    $STATS{"Len. of longest contig"} = $largestContig;
    $STATS{"Total Contig Len."} = $totalCtgLen;
    $STATS{"Average Contig Len."} = $avgCtgLen;
    $STATS{"Genome %GC Content"} = $GCcontent;
	
	
#	 if($v){
#	 
#	 $OUTPUT = "$OUTPUT"."\nRun Prefix\t$pre".
#	 "\nContig File\t$contigFile".
#	 "\nGenome Length\t$genomeLength".
#	 "\nContig count\t$contigCount".
#	 "\nContig n50\t$contig_n50".
#	 "\nContig n75\t$contig_n75".
#	 "\nContig n90\t$contig_n90".
#	 "\nLen. of longest contig\t$largestContig".
#	 "\nTotal Contig Len.\t$totalCtgLen".
#	 "\nAverage Contig Len.\t$avgCtgLen".
#	 "\nGenome %GC Content\t$GCcontent";
#	 
#	 }elsif(!$v){
#	 
#	 $OUTPUT = "$OUTPUT"."$pre".
#	 "\t$contigFile".
#	 "\t$genomeLength".
#	 "\t$contigCount".
#	 "\t$contig_n50".
#	 "\t$contig_n75".
#	 "\t$contig_n90".
#	 "\t$GCcontent".
#	 "\t$largestContig".
#	 "\t$totalCtgLen".
#	 "\t$avgCtgLen";
#	 }
#	 
	
    close(LENGTHS);    
    return $contigCount;
}

# Purpose : Calculate scaffold n50, n75, n90, largest scaffold, total and average scaffold length
# Input   : $scaffoldFile, $contigCount, $genomeLength (user-defined) 
# 
sub calculateScaffoldStats{
	
    my $scaffoldFile    = shift;
    my $contigCount     = shift;
    my $genomeLength    = shift;
    my $v               = shift;
    my ($totalScfLen,$avgScfLen,$scaffoldCount,$largestScaffold) = (0) x 4;
	
    my $tmp = seqInfo($scaffoldFile);
    my @dat = @{$tmp};
	
    # Stats come out of seqInfo in this order: (totalLength, maxScfLength, minScfLength, avgScfLength, numSequences)
    $totalScfLen     = $dat[0];
    $largestScaffold = $dat[1];
    $avgScfLen       = $dat[3];
    $scaffoldCount   = $dat[4]; 
	
	# Check if scaffold span is less than 50% of expected genome length
	#
    if($totalScfLen < ($genomeLength/2)){
		print "WARNING:  Total contig span is less than 50% of the total genome length.  Using contig span as expected genome length\n";
		printlog("WARNING:  Total contig span is less than 50% of the total genome length.  Using contig span as expected genome length\n");
		$genomeLength = $totalScfLen;
    }
	
	
	# Get nX Values
	#
    my $scaffold_n50    = getNX($scaffoldFile,50,$genomeLength);
    my $scaffold_n75    = getNX($scaffoldFile,75,$genomeLength);
    my $scaffold_n90    = getNX($scaffoldFile,90,$genomeLength);
	
	
	# Format find_nX output
	#
    chomp($scaffold_n50);
    chomp($scaffold_n75);
    chomp($scaffold_n90);
    chomp($scaffoldCount);
    chomp($largestScaffold);
    chomp($totalScfLen);
	
	# Calculate number of intrascaffold gaps
	#
    my $IntraScaffoldGaps = ($contigCount - $scaffoldCount);
    my $GapsPer5kb        = ($IntraScaffoldGaps/($genomeLength/5000));
	
	# Format Output -- if vertical option specified, print vertically.
	# if not, print horizontally
	#
	
    $STATS{"Scaffold File"} = $scaffoldFile;
    $STATS{"Scaffold count"} = $scaffoldCount;
    $STATS{"Scaffold n50"} = $scaffold_n50;
    $STATS{"Scaffold n75"} = $scaffold_n75;
    $STATS{"Scaffold n90"} = $scaffold_n90;
    $STATS{"Intra-scaffold gaps"} = $IntraScaffoldGaps;
    $STATS{"Gaps per 5kb"} = $GapsPer5kb;
	
#	=head
#	 if($v){
#	 
#	 $OUTPUT = "$OUTPUT"."\nScaffold File\t$scaffoldFile".
#	 "\nScaffold count\t$scaffoldCount".
#	 "\nScaffold n50\t$scaffold_n50".
#	 "\nScaffold n75\t$scaffold_n75".
#	 "\nScaffold n90\t$scaffold_n90".
#	 "\nIntra-scaffold gaps\t$IntraScaffoldGaps".
#	 "\nGaps per 5kb\t$GapsPer5kb";
#	 
#	 }elsif(!$v){
#	 
#	 $OUTPUT = "$OUTPUT"."\t$scaffoldFile".
#	 "\t$scaffoldCount".
#	 "\t$scaffold_n50".
#	 "\t$scaffold_n75".
#	 "\t$scaffold_n90".
#	 "\t$IntraScaffoldGaps".
#	 "\t$GapsPer5kb";
#	 }
#	 =cut
	
    close(LENGTHS);
	
    return $scaffoldCount;
}

# Purpose : Calculate read statistics on a fasta-formatted read file
# Input   : reference to a list of fragment files, vertical output flag
#
# NOTE    :  Will ONLY calculate basic stats for fragment files.
#         :  Support for mate pair stats has not been added yet.  
#
sub getReadStats{
    
    my $ref           = shift;
    my $v             = shift;
	
    my %FILES         = %{$ref};
    my %fragments     = %{$FILES{'frags'}};
    my %frmates       = %{$FILES{'frmates'}};
    my %interleaved   = %{$FILES{'interleaved'}};
    my @fragsOnly     = split(/\,/,$FILES{'frags'}{'FRAGSONLY'});
	
	
    my ($tmpA,$tmpB,$first1,$id) = "";
    my (@arrA,@arrB,@merged,@tmp) = ();
    my ($numReadsInput,$numBasesInput, $avgReadLength, $totalMates) = (0) x 4;
	
    printlog("Now analyzing reads...\n");
	
    # Count reads in the fragment-only file
    if(scalar(@fragsOnly) > 0){
		foreach my $f (@fragsOnly){
			
			if($f ne ""){
				
				printlog("Now reading: $f from fragment files ...\n");		
				chomp($f);
				
				my $tmp = seqInfo($f);
				my @dat = @{$tmp};
				
				# Count reads with sequence_info
				$numReadsInput += $dat[4];
				$numBasesInput += $dat[0];
				
				chomp($numReadsInput);
				chomp($numBasesInput);
			}
		}
    }
	
    # Go through each library, extract read files from each library, then count bases, reads and mates
    for(my $i=1;$i<=$libID;$i++){
		
		my @frags   = split(/\,/,$fragments{$i});
		my @frmates = split(/\,/,$frmates{$i});
		my @ileaved = split(/\,/,$interleaved{$i});
		my @err = ();
		
		# Put all of the read files together in order to count all bases and all reads
		@tmp   = (@frags,@frmates,@ileaved);
		
		
		# Count all fragments in libraries
		foreach my $f (@tmp){
			
			if($f ne ""){
				
				printlog("Now reading: $f from library $i ...\n");		
				chomp($f);
				
				my $tmp = seqInfo($f);
				my @dat = @{$tmp};
				
				# Count reads with sequence_info
				$numReadsInput += $dat[4];
				$numBasesInput += $dat[0];
				
				chomp($numReadsInput);
				chomp($numBasesInput);
			}
		}
		
		# Count all interleaved mates
		foreach my $readfile (@ileaved){
			if(($readfile) && (-s $readfile)){
				
				# check the beginnings of sequence file to determine fastq or fasta format
				$first1 = `head -n 1 $readfile`;
				
				# Count reads in fasta file
				if (( $first1 =~ m/^>/ )){
					
					$tmpA = `grep -ce '^>' $readfile`;
					$tmpA = ($tmpA/2);
					
					# Count reads in fastq file
				}elsif($first1 =~ m/^@/){
					
					my $id = substr($first1,0,3);
					$tmpA = `grep -ce '^$id' $readfile`;
					$tmpA = ($tmpA/2);
				}
				
				chomp($tmpA);
				$totalMates += $tmpA;
			}	
		}
		
		if (scalar(@frmates) > 0){
			
			# Assume that the number of reads in one non-interleaved mate file is the same in the other
			#    Otherwise, clc_ref_assemble will not run
			$first1 = `head -n 1 $frmates[0]`;
			
			if($first1 =~ m/^>/){
				
				$tmpA = `grep -ce '^>' $frmates[0]`;
				
			}elsif($first1 =~ m/^@/){
				
				my $id  = substr($first1,0,3);
				$tmpA = `grep -ce '^$id' $frmates[0]`;
			}
			
			chomp($tmpA);
			$totalMates += $tmpA;
		}    
    }
    
    if($numReadsInput > 0){
		
        $avgReadLength = ($numBasesInput/$numReadsInput);
        $avgReadLength = sprintf("%.2f",$avgReadLength);
		
    }
	
    # Run extra checks for empty mates
    if (!defined($totalMates)){
		$totalMates = 0;
    }
	
	# Format Output -- if vertical option specified, print vertically.
	# if not, print horizontally
	#
    $STATS{"Num. of bases in reads"} = $numBasesInput;
    $STATS{"Num. of reads used"} = $numReadsInput;
    $STATS{"Avg. read length"} = $avgReadLength;
    $STATS{"Total Mates"} = $totalMates;
	
#	=head
#	 if($v){
#	 
#	 $OUTPUT = "$OUTPUT"."\nNum. of bases in reads\t$numBasesInput".
#	 "\nNum. of reads used\t$numReadsInput".
#	 "\nAvg. read length\t$avgReadLength".
#	 "\nTotal Mates\t$totalMates\n";
#	 
#	 }elsif(!$v){
#	 
#	 $OUTPUT = "$OUTPUT"."\t$numBasesInput".
#	 "\t$numReadsInput".
#	 "\t$avgReadLength".
#	 "\t$totalMates";
#	 }
#	 =cut
	
    close(LENGTHS);
}

# Purposes : Process input files and create command strings to run clc_ref_assemble long
#          : Convert .cas assembly file to .bam assembly file for use by SAMTools
# Input    : contig file name, percent length, percent identity, hash containing options and read file names, and filename prefix
#          : (see clc_ref_assemble_long documentation for more details)
#
sub runCLCMapping{

    my $contigfile = shift;
    my $length     = shift;
    my $identity   = shift;
    my $ref        = shift;
    my $pre        = shift;

# Extract and de-reference input file names and other options from the menu hash
#
    my %READS   = %{$ref};
    
# Read files
#
    my %frags    = %{$READS{'frags'}};
    my %frmates  = %{$READS{'frmates'}};  # Separate forward and reverse mate files
    my %sffflx   = %{$READS{'sffflx'}};
    my %sfftitan = %{$READS{'sfftitan'}};
    my %ileaved  = %{$READS{'interleaved'}};
    my @fragsOnly  = split(/\,/,$frags{'FRAGSONLY'});
    
# Paired End/Mate Options 
#
    my %distmode = %{$READS{'distmode'}};
    my %orient   = %{$READS{'orient'}};
    my %insert   = %{$READS{'insert'}};
    
# Command strings
#
    my $frag_cmd = "";
    my $mate_cmd = "";
    
    printlog("Mapping reads to contigs...\n ");
    printlog("Aligning reads using contig file: ".$contigfile."\n");
    
# Go through each insert library, Append filenames and options to command strings that get passed to the clc_ref_assebler
#
    for(my $i=1;$i<=$libID;$i++){

	#VERY IMPORTANT: Clear the fragment command for each library!
	$frag_cmd = "";
	
	my (@frags,@frmates,@insert,@fragsOnly) = ();
	my ($orient,$distmode,$ileaved,$tmp) = "";

	@frags     = split(/\,/,$frags{$i});
	@fragsOnly = split(/\,/,$frags{'FRAGSONLY'});

	if((scalar(@frags)) > 0){
	    
	    # Check fragment files
	    foreach my $file (@frags){
#=head
#                This isn't necessary, since checkFiles already found out whether or not all input
#		were present and not empty
#
#		if(!(-e $file)){die "\nFragment file $file does not exist!\n\n";}
#		elsif(-z $file){die "\nFragment file $file is empty!\n\n";}
#=cut
		$frag_cmd = $frag_cmd." ".$file;
	    }
	    
	}else{printlog("No fragment files given.\n");}


# If we have mates, append the mate options to the command string (plus clc flags)
#   Do not run mate mapping if insert lengths are not given!
#
	@frmates  = split(/\,/,$frmates{$i});
	@insert   = split(/\,/,$insert{$i});
	$ileaved  = $ileaved{$i};
	$orient   = $orient{$i};
	$distmode = $distmode{$i};

	if(!$orient){
	    printlog("Mate pair orientation not specified.  Using ss as default.\n");
	    $orient =  $DEFAULT_DIST_MODE;
	}
	if(!$distmode){
	    printlog("Distance mode not specified.  Using ee as default.\n");
	    $distmode = $DEFAULT_DIRECTION;
	}
	
	if(scalar(@frmates) > 0){
	    foreach my $file (@frmates){

#		if(!(-e $file)){die "\nMate file $file does not exist!\n\n";}
#		elsif(-z $file){die "\nMate file $file is empty!\n\n";}
	    }
	    $tmp = " ".join(' ',@frmates);
	}
	
	if($ileaved){

#	    if(-z $ileaved){die "\nMate file $ileaved is empty!\n\n";}

	    $tmp = " ".$ileaved;
	}
	$mate_cmd = $mate_cmd.$tmp if(defined($tmp));
    }

    if( scalar(@fragsOnly)> 0){
	$frag_cmd = $frag_cmd." ".join(" ",@fragsOnly);
    }

# Finally... run the clc mapper
#
#    system("$clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd");
	 printl("++ RUN CMD : $clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd 2>&1\n");
     my $tolog = `$clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd 2>&1`;
    
    
# Write status to logfile
#
    printlog("CLC_COMMAND = $clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd\n");
    printlog("\n### CLC REF ASSEMBLE OUTPUT  ###");

    # Check if clc_ref_assemble_long had any problems opening any of the sequence files
    if($tolog =~ m/(P|p)roblem/){
	my @errs = split(/\s+/,$tolog);
	my $err  = pop(@errs);
	printlog("Problem opening file $err.  Please verify file format and location.\n");
	die "Problem opening file $err.  Please verify file format and location.\n";
    }

    printlog($tolog);
    printlog("\n################################\n");
    printlog("\nMapping Complete.\n");

    
# Check if assembly was completed by reading the .cas file.  Will be empty if nothing assembled.
# (Written for debug purposes)
#
    if(-e "$pre.cas"){
	if(-z "$pre.cas"){ 
	    die "\nRead mapping incomplete.  Please check the CLC_COMMAND in the logfile.\n\n";
	}
    }else{ 
	die "\nRead mapping incomplete.  Please check the CLC_COMMAND in the logfile.\n\n";
    }
    
# convert CAS assembly to BAM format
#
    printlog("Converting $pre.cas to $pre.bam...\n");
    printl("++ RUN CMD : $castosam -a $pre.cas -o $pre.bam 2>> samtools.err\n");
    my $err = `$castosam -a $pre.cas -o $pre.bam 2>> samtools.err`;
    

# Clean-up big output files
#
    if(-e "./$prefix.frags.fasta"){system "/bin/rm ./$prefix.frags.fasta";}
    if(-e "./$prefix.mates.fasta"){system "/bin/rm ./$prefix.mates.fasta";}
    if(-e "./$pre.cas"){
	if($keep){
	    printlog(".cas file stored in working directory\n");
	}else{
	    system "/bin/rm ./$pre.cas";
	    printlog(".cas assembly file removed!\n");
	}
    }

}
#sub runCLCMapping{
#	
#    my $contigfile = shift;
#    my $length     = shift;
#    my $identity   = shift;
#    my $ref        = shift;
#    my $pre        = shift;
#	
#	# Extract and de-reference input file names and other options from the menu hash
#	#
#    my %READS   = %{$ref};
#    
#	# Read files
#	#
#    my %frags    = %{$READS{'frags'}};
#    my %frmates  = %{$READS{'frmates'}};  # Separate forward and reverse mate files
#    my %sffflx   = %{$READS{'sffflx'}};
#    my %sfftitan = %{$READS{'sfftitan'}};
#    my %ileaved  = %{$READS{'interleaved'}};
#    my @fragsOnly  = split(/\,/,$frags{'FRAGSONLY'});
#    
#	# Paired End/Mate Options 
#	#
#    my %distmode = %{$READS{'distmode'}};
#    my %orient   = %{$READS{'orient'}};
#    my %insert   = %{$READS{'insert'}};
#    
#	# Command strings
#	#
#    my $frag_cmd = "";
#    my $mate_cmd = "";
#    
#    printlog("Mapping reads to contigs...\n ");
#    printlog("Aligning reads using contig file: ".$contigfile."\n");
#    
#	# Go through each insert library, Append filenames and options to command strings that get passed to the clc_ref_assebler
#	#
#    for(my $i=1;$i<=$libID;$i++){
#		
#		my (@frags,@frmates,@insert,@fragsOnly) = ();
#		my ($orient,$distmode,$ileaved,$tmp) = "";
#		
#		@frags     = split(/\,/,$frags{$i});
#		@fragsOnly = split(/\,/,$frags{'FRAGSONLY'});
#		
#		if((scalar(@frags)) > 0){
#			
#			# Check fragment files
#			foreach my $file (@frags){
#
#				$frag_cmd = $frag_cmd." ".$file;
#				printlog("Fragment file $file added to fragment portion of command.\n");
#			}
#			
#		}else{printlog("No fragment files given.\n");}
#		
#		
#		# If we have mates, append the mate options to the command string (plus clc flags)
#		#   Do not run mate mapping if insert lengths are not given!
#		#
#		@frmates  = split(/\,/,$frmates{$i});
#		@insert   = split(/\,/,$insert{$i});
#		$ileaved  = $ileaved{$i};
#		$orient   = $orient{$i};
#		$distmode = $distmode{$i};
#		
#		if(!$orient){
#			printlog("Mate pair orientation not specified.  Using fb as default.\n");
#			$orient =  $DEFAULT_DIST_MODE;
#		}
#		if(!$distmode){
#			printlog("Distance mode not specified.  Using ss as default.\n");
#			$distmode = $DEFAULT_DIRECTION;
#		}
#		
#		if(scalar(@frmates) > 0){
#			foreach my $file (@frmates){
#				
#				#		if(!(-e $file)){die "\nMate file $file does not exist!\n\n";}
#				#		elsif(-z $file){die "\nMate file $file is empty!\n\n";}
#			}
#			$tmp = " -p $orient $distmode ".join(" ",@insert)." -i ".join(' ',@frmates);
#		}
#		
#		if($ileaved){
#			
#			#	    if(-z $ileaved){die "\nMate file $ileaved is empty!\n\n";}
#			
#			$tmp = " -p $orient $distmode ".join(" ",@insert)." ".$ileaved;
#		}
#		$mate_cmd = $mate_cmd.$tmp if(defined($tmp));
#    }
#	
#    if( scalar(@fragsOnly)> 0){
#		$frag_cmd = $frag_cmd." ".join(" ",@fragsOnly);
#    }
#	
#	# Finally... run the clc mapper
#	#
#	#    runsys("$clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd");
#	printl("++ RUN CMD : $clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd 2>&1\n");
#	my $tolog = `$clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd 2>&1`;
#    
#    
#	# Write status to logfile
#	#
#    printlog("CLC_COMMAND = $clc_ref_assemble_long -o $pre.cas -d $contigfile -l $length -s $identity -r random -q $frag_cmd $mate_cmd\n");
#    printlog("\n### CLC REF ASSEMBLE OUTPUT  ###");
#	
#    # Check if clc_ref_assemble_long had any problems opening any of the sequence files
#    if($tolog =~ m/(P|p)roblem/){
#		my @errs = split(/\s+/,$tolog);
#		my $err  = pop(@errs);
#		printlog("Problem opening file $err.  Please verify file format and location.\n");
#		die "Problem opening file $err.  Please verify file format and location.\n";
#    }
#	
#    printlog($tolog);
#    printlog("\n################################\n");
#    printlog("\nMapping Complete.\n");
#	
#    
#	# Check if assembly was completed by reading the .cas file.  Will be empty if nothing assembled.
#	# (Written for debug purposes)
#	#
#    if(-e "$pre.cas"){
#		if(-z "$pre.cas"){ 
#			die "\nRead mapping incomplete.  Please check the CLC_COMMAND in the logfile.\n\n";
#		}
#    }else{ 
#		die "\nRead mapping incomplete.  Please check the CLC_COMMAND in the logfile.\n\n";
#    }
#    
#    
#    
#    
#	
#	# Clean-up big output files
#	#
#	#    if(-e "./$prefix.frags.fasta"){system "/bin/rm ./$prefix.frags.fasta";}
#	#    if(-e "./$prefix.mates.fasta"){system "/bin/rm ./$prefix.mates.fasta";}
#    if(-e "./$pre.cas"){
#		if($keep){
#			printlog(".cas file stored in working directory\n");
#		}else{
#			#	    system "/bin/rm ./$pre.cas";
#			printlog(".cas assembly file removed!\n");
#		}
#    }
#	
#}

##### 7/18 - BAD GREP METHOD
# Purposes : Calculate average coverage, individual contig coverage and % of genome with coverage > 5X from the assembly BAM file
# Input    : contig file name, filename prefix, vertical output flag 
#
#sub coverageStats{
#	my $contigFile = shift;
#    my $pre        = shift;
#    my $v          = shift;
#    my $covCut     = shift;
#	
#    my ($basecount,$gt5X,$percGt5X,$totalCov,$avgCov) = (0) x 5;
#    my $bamfile    = "$pre.bam";
#    my %LENG;
#	my @contigNamesOrdered;
#	
#    ## QUICK FIX! ##
#                    	printl("++ RUN CMD : $sequence_info -l $contigFile\n");
#    my @out = `$sequence_info -l $contigFile`;    
#    @out = @out[13..($#out)];
#    foreach my $l(@out) {
#		$l =~ s/^\s+//;
#		chomp($l);
#		my @dat = split(/\s+/,$l);
#		my $len = $dat[$#dat];
#		my $id  = $dat[0];
#		
#		if((defined($id)) && (defined($len))) {
#			$LENG{$id} = $len;
#			push(@contigNamesOrdered,$id);
#		}
#    }
#	
#	# Execute assembly_info to acquire each contig's average coverage.
#	printl("++ RUN CMD : $assembly_info $pre.cas 2>> assembly_info.error\n");
#	my @ai_out = `$assembly_info $pre.cas 2>> assembly_info.error`;    
#    @ai_out = @ai_out[35..($#ai_out)-1];
#    foreach my $l (@ai_out) {
#		$l =~ s/^\s+//;
#		chomp($l);
#		my @dat = split(/\s+/,$l);
#		my $id  = @contigNamesOrdered[int($dat[0])-1];
#		my $cov = $dat[-1];
#		my $len = 0;
#		
#		if((defined($id)) && (defined($cov))) {
#			printlog("\t $id\_AVERAGE CONTIG COVERAGE = $cov\n");
#			$len = $LENG{$id} if (defined($LENG{$id}));
#			$totalCov += ($cov * $len);
#		}
#    }
#	
#	# Execute assembly_table to acquire percentage of genome with coverage at covCut or higher.
#	printlog("Acquiring assembly table of cas file\n\n");
#    runsys("$assembly_table $pre.cas > assembly_table.txt 2>> assembly_table.error");
#	
#	open(ALIGNMENTS,"<assembly_table.txt") or die "Problem opening assembly_table file.\n $!";
#	open(ALIGNMENTS_BY_BASE,">assembly_table_by_base.txt") or die "Problem opening assembly_table_by_base file.\n $!";
#    
#	while (my $line = <ALIGNMENTS>) {
#		my @fields = split(/\s+/,$line);
#		next if (int($fields[4]) == -1);
#		
#		my $contig = @contigNamesOrdered[int($fields[4])-1];
#		my $refStart = int($fields[5]);
#		my $refEnd = int($fields[6]);
#		
#		for (my $count = $refStart; $count <= $refEnd; $count++) {
#			print ALIGNMENTS_BY_BASE "$contig,$count\n";
#		}
#	}
#    
#	close(ALIGNMENTS);
#	close(ALIGNMENTS_BY_BASE);
#	
#	printl("\nGrep CTG Cov - Start\n");
#	supply_time();
#	
#    foreach my $contig (@contigNamesOrdered) {
#		my $contigLength = $LENG{$contig};
#		
#		for (my $baseCoord = 1; $baseCoord <= $contigLength; $baseCoord++) {
#			
#			printl("++ RUN CMD : grep -c \"$contig,$baseCoord\" assembly_table_by_base.txt\n");
#			my $baseCov = `grep -c "$contig,$baseCoord" assembly_table_by_base.txt`;
#			runsys("echo -e \"$contig, $baseCoord, $baseCov\" >> listings.txt");
#			$gt5X++ if (int($baseCov) >= $covCut);
#		}
#		
#		$basecount += $contigLength;
#	}
#	printl("\nGrep CTG Cov - End\n");
#	supply_time();
#	
#    if ($basecount > 0) {
#		$avgCov = ($totalCov/$basecount);    
#		
#		if($gt5X < $basecount) {
#			$percGt5X = (($gt5X/$basecount)*100);
#		}
#		else {
#			print "ERROR: Number of bases with > 5x coverage greater than genome length.\nPlease check contig and read files.\n";
#			printlog("ERROR: Number of bases with > 5x coverage greater than genome length.\nPlease check contig and read files.\n");
#		}
#		
#		$avgCov   = sprintf("%.2f",$avgCov);
#		$percGt5X = sprintf("%.2f",$percGt5X);
#    }
#	else {
#		print "ERROR: No bases found in coverage file\n";
#		printlog("ERROR: No bases found in coverage file\n");
#    }
#	
#	# Format Output -- if vertical option specified, print vertically.
#	# if not, print horizontally
#	#
#    $STATS{"Average genome coverage"} = $avgCov;
#    $STATS{"Percent coverage"} = $percGt5X;
#	
#}


sub coverageStats{

    my $contigFile = shift;
    my $pre        = shift;
    my $v          = shift;
    my $covCut     = shift;

    my ($basecount,$gt5X,$percGt5X,$totalCov,$avgCov) = (0) x 5;
    my $bamfile    = "$pre.bam";
    my %CTG;
    my %LENG;

    ## QUICK FIX! ##
    my @out = `$sequence_info -l $contigFile`;    
    @out = @out[13..($#out)];
    foreach my $l(@out){

	$l =~ s/^\s+//;
	chomp($l);
	my @dat = split(/\s+/,$l);
	my $len = $dat[$#dat];
	my $id  = $dat[0];

	if((defined($id)) && (defined($len))){$LENG{$id} = $len;}
    }

    printlog("\nNow Sorting: $bamfile to $bamfile.sorted \n");
    runsys("$SAMTOOLS_DIR/samtools sort $bamfile $bamfile.sorted 2>> samtools.err");


# Add the .sam file extension because samtools sort adds the bam suffix AFTER the "sorted" suffix
#
    $bamfile = $bamfile."\.sorted";
    $bamfile = $bamfile."\.bam";
    printlog("Now piling: $bamfile \n\n");
    runsys("$SAMTOOLS_DIR/samtools mpileup -f $contigFile $bamfile > ./coverage.txt 2>> samtools.err");

# Create a hash with contig ID as the key 
# and a list of coverage values for that contig as the value
#
    open(BASES,"<coverage.txt") or die "Problem opening coverage file.\n $!";
    while(my $line = <BASES>){

	my @fields = split(/\s+/,$line);
	my @covg;

	if(exists($CTG{$fields[0]})){
	    unshift @{$CTG{$fields[0]}},$fields[-3];
	}else{
	    unshift @covg,$fields[-3];
	    $CTG{$fields[0]} = \@covg;
	}
	if($fields[-3] >= $covCut){$gt5X ++;}
    }

# Cleanup
#
    close(BASES);
#    if(-e "./coverage.txt"){system "/bin/rm ./coverage.txt"}

# Calculate contig coverage as a "mean of means"
# This means getting the average coverage for each contig
# then calculate the overall average coveage by averaging the average contig coverage
#
    foreach my $contig(keys%CTG){

	my ($tmp,$pos,$contig_cov,$len) = 0;

	## SPECIAL CASE: A contig may have no reads mapped to it
	if(defined($LENG{$contig})){$len = $LENG{$contig};}else{$len=0;}

	$basecount += $len;

	foreach my $cov(@{$CTG{$contig}}){    

	    $tmp += $cov;
	    $pos ++;
	}
	
	$contig_cov = ($tmp/$pos);	
	$contig_cov = sprintf("%.2f",$contig_cov);

	printlog("\t $contig\_AVERAGE CONTIG COVERAGE = $contig_cov\n");
	printlog("\t $contig\_CONTIG LENGTH = $pos\n");
	
	$totalCov += ($contig_cov * $len);
    }

# Divide the total coverage per contig by the number of contigs
# NOTE: using basecount here may not be correct.  Use contig lengths derived from sizefasta or SAMtools?
#
    if($basecount > 0){

	$avgCov   = ($totalCov/$basecount);    

	if($gt5X < $basecount){
	    $percGt5X = (($gt5X/$basecount)*100);
	}else{
	    print "ERROR: Number of bases with > 5x coverage greater than genome length.\nPlease check contig and read files.\n";
	    printlog("ERROR: Number of bases with > 5x coverage greater than genome length.\nPlease check contig and read files.\n");
	}

	$avgCov   = sprintf("%.2f",$avgCov);
	$percGt5X = sprintf("%.2f",$percGt5X);

    }else{
	
	print "ERROR: No bases found in coverage file\n";
	printlog("ERROR: No bases found in coverage file\n");
    }

# Format Output -- if vertical option specified, print vertically.
# if not, print horizontally
#
    $STATS{"Average genome coverage"} = $avgCov;
    $STATS{"Percent coverage"} = $percGt5X;

#=head
#    if($v){
#
#	$OUTPUT = "$OUTPUT"."\nAverage genome coverage\t$avgCov".
#	    "\nPercent of genome with >= $covCut X coverage\t$percGt5X";
#
#    }elsif(!$v){
#
#	$OUTPUT = "$OUTPUT"."\t$avgCov".
#	    "\t$percGt5X";
#    }
#=cut

# Cleanup
#
    # Alwasys remove original bam file
    if(-e "./$pre.bam"){
       	    printlog("Removed $pre.bam!\n");
#	    system "/bin/rm ./$pre.bam";
    }

    # Keep the sorted copy of the bam file, if needed
    if(-e "./$pre.bam.sorted.bam"){
	if($keep){
	    printlog("$pre.bam.sorted.bam stored in current working directory\n");
	}else{
	    printlog("Removed $pre.bam.sorted.bam!\n");
#	    system "/bin/rm ./$pre.bam.sorted.bam";
	}
    } 

    # Keep the fasta index, as it may be needed for other 
    # scripts to parse the bam/sam file
    if(-e "./$pre.fai"){
	if($keep){
	    printlog("$pre.fai stored in current working directory\n");
	}else{
	    printlog("Removed $pre.fai!\n");
#	    system "/bin/rm ./$pre.fai";
	}
    }
}


sub seqInfo{
    
    my $file   = shift;
    my $seqInf = "/usr/local/packages/clc-ngs-cell/sequence_info";
    my ($tot,$max,$min,$avg,$seq) = (0) x 5;
    my @dat;
    
    # run sequence_info, capture STDOUT and STDERR
    printl("++ RUN CMD : $seqInf $file 2>&1\n");
    my @infos = `$seqInf $file 2>&1`;
    
    
    foreach my $l(@infos){
		
		# sequence_info prints the word "Problem" or "problem" when there is 
		# an issue with the format of the files given to it
		if($l =~ m/(P|p)roblem/){ 
			printlog("Problem with format of fasta file in $file\t");
			die "Problem with format of fasta file in $file\t";  
		}
		
		my @t = split(/\s+/,$l);
		if(defined($t[1])){
			if( $t[1] =~ m/Total/   ){ $tot = $t[2]; }
			if( $t[1] =~ m/Maximum/ ){ $max = $t[2]; }
			if( $t[1] =~ m/Minimum/ ){ $min = $t[2]; }
			if( $t[1] =~ m/Average/ ){ $avg = $t[2]; }
			if( $t[0] =~ m/Number/  ){ $seq = $t[3]; }
		}
    }
    
    chomp($tot);
    chomp($max);
    chomp($min);
    chomp($avg);
    chomp($seq);
    
    @dat = ($tot,$max,$min,$avg,$seq);
    
    return \@dat;
}
			
################################################################################################################################################################################################################################################
sub ZeroXAnalysis {		
# ZeroXAnalysis($contigFile,$prefix,$zero_x_ref);
	my $contigFile   = shift;
	my $prefix   = shift;
	my $zero_x_ref   = shift;
	
	# clc_ref_assemble_long -o INKIN.orig_best_spades.newb_shred.cas -d /local/netapp_scratch/CORE/jmccorri/EUK/REFERENCE/t_asahii_ref.fasta -q ./INKIN.orig_best_spades.newb_shred.fasta
	runsys("ln -s /usr/local/packages/clc-ngs-cell/license.properties .");
	runsys("/usr/local/devel/BCIS/assembly/tools/contig_quik_shred.pl $contigFile 7999 49 ./zeroX.shred_ctg.fasta");
	runsys("/usr/local/packages/clc-ngs-cell/clc_ref_assemble_long -o zeroX.cas -d $zero_x_ref -q ./zeroX.shred_ctg.fasta");
	# assembly_info -c 01_PE_only.newb_shred.cas | grep 'covered 0' | head -n 1
	printl("/usr/local/packages/clc-ngs-cell/assembly_info -c zeroX.cas | grep \'covered 0\' | head -n 1 | awk \'{print \$5}\'\n");
	my $return = `/usr/local/packages/clc-ngs-cell/assembly_info -c zeroX.cas | grep \'covered 0\' | head -n 1 | awk \'{print \$5}\'`;
	
	my $tmp = seqInfo($zero_x_ref);
	my @dat = @{$tmp};
	my $genomespan = $dat[0];
	
	my $perc = sprintf( "%.2f", (($return/$genomespan)*100) ) ;
	
	return $perc;
}
			
################################################################################################################################################################################################################################################
sub supply_time {
	
	#requires declared globals @prev_time and $time_bool=0
	
	my @weekday = ("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat");
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
	my @cur_time = ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
	$year = $year + 1900;
	$mon += 1;
	
	
	if ($time_bool == 0){
		#first time track
	}
	else{
		#compare
		
		#print "prev_time 0 (A) = 0 $prev_time[0] 1 $prev_time[1] 2 $prev_time[2] 3 $prev_time[3] 4 $prev_time[4] 5 $prev_time[5]\n";
		
		my $prev_sec = $prev_time[0];
		my $prev_min = $prev_time[1];
		my $prev_hour = $prev_time[2];
		my $prev_mday = $prev_time[3];
		my $prev_mon = $prev_time[4] + 1;
		my $prev_year = $prev_time[5] + 1900;
		
		my $seconds=0;
		
		#print "year : $year / $prev_year\n";
		#print "mon : $mon / $prev_mon\n";
		#print "mday : $mday / $prev_mday\n";
		
		if (($year == $prev_year) && ($mon == $prev_mon) && ($mday == $mday)){
			#SAME DAY
				my $cur_daytime = (($hour*3600)+($min*60)+$sec);
				my $prev_daytime = (($prev_hour*3600)+($prev_min*60)+$prev_sec);
				my $diff = $cur_daytime - $prev_daytime;
				printl("+++ TIME : $mday/$mon/$year $hour:$min:$sec $weekday[$wday]\n");
				#printl("( $cur_daytime - $prev_daytime )\n");	
				printl("+++ TIME EXPIRED : $diff seconds\n");
		}else{
			#DIFFERENT DAY
			if ($mday > $prev_mday){
				#WITHIN A MONTH
				my $day_diff = (($mday-$prev_mday) * 86400);
				my $cur_daytime = (($hour*3600)+($min*60)+$sec);
				my $prev_daytime = (($prev_hour*3600)+($prev_min*60)+$prev_sec);
				my $diff = $day_diff + $cur_daytime + (86400 - $prev_daytime);
				printl("+++ TIME : $mday/$mon/$year $hour:$min:$sec $weekday[$wday]\n");	
				printl("+++ TIME EXPIRED : $diff seconds\n");
			}else{
				#OUTSIDE A MONTH
				printl("+++ TIME : $mday/$mon/$year $hour:$min:$sec $weekday[$wday]\n");	
				printl("+++ TIME EXPIRED : CROSSES MONTH BOUNDS!  SEE TIMESTAMP!\n");
			}
		}
	}
	
	
	$time_bool = 1;
	@prev_time = @cur_time;
	
		#print "prev_time 0 (B) = $prev_time[0]\n";
}

			
			
			################################################################################################################################################################################################################################################
			# runsys : 
			################################################################################################################################################################################################################################################
			sub runsys {
				my $cmd = shift;
				print "++ RUN CMD : $cmd\n";
				print LOG "++ RUN CMD : $cmd\n";
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
			

		
### END OF SCRIPT ###
