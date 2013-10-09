#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2004 and later by Bastien Chevreux
#
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the 
# Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
# 
#

#set min 3
#set max 7
#for {set i 0} {$i< 100} {incr i} {
#    set nummut [expr {$min+int(rand()*($max-$min+1))}]
#    puts $nummut
#}
#exit

namespace eval fasta2frag {
  variable opts
}

proc string_reverse str {
    set res {}
    set i [string length $str]
    while {$i > 0} {append res [string index $str [incr i -1]]}
    set res
}

proc list_reverse list {
    set _temp {}
    for	{set i [ expr [ llength $list ] - 1 ] } {$i >= 0} {incr i -1}  	{
	lappend _temp [ lindex $list $i ]
    }
    return  $_temp
}


proc fasta2frag::dnamut_basechange {dnaseq} {
    set lenseq [string length $dnaseq]
    set cpos [expr {int(rand()*$lenseq)}]
    set actbase [string index $dnaseq $cpos]
    set tmp [expr {int(rand()*3)}]

    set i 0
    foreach base {a c g t} {
	if {$base != $actbase} {
	    if {$i == $tmp} {
		set newbase [string toupper $base]
	    }
	    incr i
	}
    }
    
    puts "oldseq: $dnaseq"
    set dnaseq [string replace $dnaseq $cpos $cpos $newbase]
    puts "newseq: $dnaseq"
    return $dnaseq
}

proc fasta2frag::dnamut_insert {dnaseq} {
    set lenseq [string length $dnaseq]
    set inspos [expr {int(rand()*$lenseq)}]
    set tmp [expr {int(rand()*4)}]
    set inswhat "[string index $dnaseq $inspos][lindex {a c g t} $tmp]"
    puts "insert: $insbase"
    set dnaseq [string replace $dnaseq $inspos $inspos $inswhat]
    return $dnaseq
}

proc fasta2frag::dnamut_delete {dnaseq} {
    set lenseq [string length $dnaseq]
    set delpos [expr {int(rand()*$lenseq)}]
    set dnaseq [string replace $dnaseq $delpos $delpos]
    return $dnaseq
}

proc fasta2frag::dnamut {dnaseq min max} {
    set nummut [expr {$min+int(rand()*($max-$min+1))}]

    for {set i 1} {$i <= $nummut} {incr i} {
	puts "mut $i"
	#set dnaseq [dnamut_basechange $dnaseq]
	
	#set muttype [expr {int(rand()*3)}]
	#switch -exact -- $muttype {
	#    0 { set $dnaseq [dnamut_basechange $dnaseq] }
	#    1 { set $dnaseq [dnamut_insert $dnaseq] }
	#    2 { set $dnaseq [dnamut_delete $dnaseq] }
	#}
    }

    return $dnaseq
}

proc fasta2frag::dumpfasta {name descline dnaseq fout qualseq qout} {
    variable opts

    if {$opts(-minmut) || $opts(-maxmut)} {
	set dnaseq [dnamut $dnaseq $opts(-minmut) $opts(-maxmut)]
    }
	
    puts $fout ">$name $descline\n$dnaseq"
    if {[llength $qualseq]} {
	puts $qout ">$name\n$qualseq"
    }
}

proc fasta2frag::getStrobeData {dnaseq quals from len} {
    variable opts

    set retseq ""
    set retquals {}

    set dnalen [string length $dnaseq]

    set thisstrobe $opts(-strobeon)
    set totalbases 0
    set ii 0

    for {} {$totalbases < $len} {incr ii} {
	incr totalbases
	set dnapos [expr {($from+$ii)}]

	if {$dnapos < $dnalen || $opts(-c)} {
	    set realindex [expr {$dnapos%$dnalen}]
	    append retseq [string index $dnaseq $realindex ]
	    if {[llength $quals]>0} {
		lappend retquals [lindex $quals $realindex]
	    } else {
		lappend retquals $opts(-q)
	    }
	    incr thisstrobe -1
	    if {$thisstrobe ==0 && $totalbases < $len} {
		if {1} {
		    # fixed number of 'N', +/- 20% real advance
		    # this is how I expect PacBio data to be
		    for {set i 0} {$i < $opts(-strobeoff)} {incr i} {
			append retseq n
			lappend retquals 0
		    }
		    # variance of +/- 20%
		    set iiincr [expr {int($opts(-strobeoff)+(rand() * $opts(-strobeoff) * .4)-($opts(-strobeoff)*.2))}]
		    incr ii $iiincr
		    set thisstrobe $opts(-strobeon)
		} else {
		    # +/- 20% number of 'N', fixed real advance
		    # for testing only, this should never be a real world case!
		    
		    set fto [expr {int($opts(-strobeoff)+(rand() * $opts(-strobeoff) * .4)-($opts(-strobeoff)*.2))}]
		    for {set i 0} {$i < $fto} {incr i} {
			append retseq n
			lappend retquals 0
		    }
		    
		    incr ii $opts(-strobeoff)
		    set thisstrobe $opts(-strobeon)
		}
	    }
	}
    }

    return [list $retseq $retquals]
}

#proc fasta2frag::getdnafromto {dnaseq from to} {
#    variable opts
#
#    set retseq ""
#
#    set dnalen [string length $dnaseq]
#
#    set totalbases 0
#    set len [expr {$from - $to}]
#    set ii 0
#
#    for {} {$totalbases < $len} {incr ii} {
#	set realindex [expr {($from+$ii)%$dnalen}]
#
#	append retseq [string index $dnaseq $realindex ]
#	if {[llength $quals]>0} {
#	    lappend retquals [lindex $quals $realindex]
#	} else {
#	    lappend retquals $opts(-q)
#        }
#    }
#}

proc fasta2frag::getqualfromto {dnaseq quals from to} {
    variable opts
    
    set qualslen [llength $quals]

    set qualseq {}
    if {$qualslen > 0} {
	for {set i $from} {$i < $to} {incr i} {
	    set thisqual -1
	    catch {
		set thisqual [expr {int([lindex $quals ($i%$qualslen) ] / $opts(-qualdivisor))}]
	    }
	    
	    # noqual?
	    if {$thisqual < 0} {
		set thisqual $opts(-minqual)
	    }
	    
	    if {$thisqual < $opts(-minqual)} {
		set thisqual $opts(-minqual)
	    }
	    lappend qualseq $thisqual
	}
    } else {
	for {set i $from} {$i < $to} {incr i} {
	    lappend qualseq $opts(-q)
        }
    }

    if {[llength $qualseq] > [string length $dnaseq]} {
	set qualseq [lrange $qualseq 1 [string length $dnaseq]]
    }

    return $qualseq
}

proc fasta2frag::reverse_seq_qual {s q} {
    upvar $s dnaseq
    upvar $q qualseq
    set dnaseq [string map {a t A T c g C G g c G C t a T A } [string_reverse $dnaseq]]
    set qualseq [list_reverse $qualseq]
}

proc fasta2frag::fragASeq {nonamedata name fout quals qout mpout mbout} {
    variable opts
    #set fid [open $opts(-infile) r]
    #set data [read $fid]
    #close $fid
    
    puts "fragging $name"
    
    #if {[regsub -all -line {^>(.*?)$} $data "" nonamedata] >1} {
    #	puts stderr "More than one sequence in fasta file, aborting."
    #	exit 500
    #}
    
    regsub -all {\n} $nonamedata "" data
    unset nonamedata
    
    # when circularising is on, append first $opts(-l)-1 bases
    #  of sequence to end of sequence (and also quals if present)
    # if strobing, also add $opts(-l)-$opts(-strobeoff) bases
    if {$opts(-c) > 0} {
	puts -nonewline "Circularising ..."
	set endpos [expr {$opts(-l)-2}]

	if {$opts(-s) > 0} {
	    puts -nonewline " with strobe settings ..."
	    set endpos [expr {$endpos + ($opts(-l)-2 - $opts(-strobeoff))}]
	}

	#puts [string range $data 0 $endpos]
	append data [string range $data 0 $endpos]
	if {[llength $quals] > 0} {
	    set quals [concat $quals [lrange $quals 0 $endpos]]
	}

	puts " done."
    }

    set datalen [string length $data]

    set from 0
    set to $opts(-l)
    if {$to > $datalen} {
	set to $datalen
    }

    switch -exact -- $opts(-pairednaming) {
	454 {
	    set fwdname "_pe.f"
	    set revname "_pe.r"
	}
	solexa {
	    set fwdname "_pe/1"
	    set revname "_pe/2"
	}
	tigr {
	    set fwdname "TF"
	    set revname "TR"
	}
	default {
	    set fwdname "_pe.p1"
	    set revname "_pe.q1"
	}
    }

    set seqid 1

    # of we reach the end of the sequence, this will be set to 1
    set lastflag 0

    proc putz {bla} {
	puts $bla
    }

    for { } {$from < $datalen} { incr from $opts(-i); incr to $opts(-i); incr seqid} {
	
	if {$lastflag == 0} {
	    #puts "$from $to"
	    
	    if {$to >= $datalen} {
		set lastflag 1
	    }
	    
	    if {$opts(-p) > 0 } {
		set pairnamef ${name}_${seqid}${fwdname}
		set pairnamer ${name}_${seqid}${revname}
		if { $opts(-s) > 0 } {
		    set to2 [expr {$from + $opts(-insert_size)}]

		    set fractionoff [expr {1.0/$opts(-strobeon)*$opts(-strobeoff)}]
		    set mult [expr {1.0+$fractionoff}]
		    set from2 [expr {int($to2 - ($opts(-l)*$mult)) }]

		    if {$to2 <= $datalen} {

			#putz "F: $from $opts(-l)"
			set bla [getStrobeData $data $quals $from $opts(-l)]
			set dnaseq [lindex $bla 0]
			set qualseq [lindex $bla 1]
			
			set descline "[expr {$opts(-insert_size)-$opts(-insert_stdev)}] $opts(-insert_size) [expr {$opts(-insert_size)+$opts(-insert_stdev)}] 1 [string length $dnaseq]"

			dumpfasta $pairnamef $descline $dnaseq $fout $qualseq $qout
			
			#putz "R: $from2 $opts(-l)"
			set bla [getStrobeData $data $quals $from2 $opts(-l)]
			set dnaseq [lindex $bla 0]
			set qualseq [lindex $bla 1]
			
			reverse_seq_qual dnaseq qualseq
			dumpfasta $pairnamer $descline $dnaseq $fout $qualseq $qout

			puts $mpout "$pairnamef\t$pairnamer"
			puts $mbout "$pairnamef\t$pairnamer\tsimplepair"
		    }
		} else {
		    set to2 [expr {$from + $opts(-insert_size)}]
		    set from2 [expr {$to2 - $opts(-l)}]
		    if {$to2 <= $datalen} {
			
			set dnaseq [string range $data $from [expr {$to-1}]]
			set qualseq [getqualfromto $dnaseq $quals $from $to]

			set descline "[expr {$opts(-insert_size)-$opts(-insert_stdev)}] $opts(-insert_size) [expr {$opts(-insert_size)+$opts(-insert_stdev)}] 1 [string length $dnaseq]"
			dumpfasta $pairnamef $descline $dnaseq $fout $qualseq $qout
			
			set dnaseq [string range $data $from2 [expr {$to2-1}]]
			set qualseq [getqualfromto $dnaseq $quals $from2 $to2]
			reverse_seq_qual dnaseq qualseq
			dumpfasta $pairnamer $descline $dnaseq $fout $qualseq $qout

			puts $mpout "$pairnamef\t$pairnamer"
			puts $mbout "$pairnamef\t$pairnamer\tsimplepair"
		    }
		}
	    } else {
		if { $opts(-s) > 0 } {
		    set bla [getStrobeData $data $quals $from $opts(-l)]
		    set dnaseq [lindex $bla 0]
		    set qualseq [lindex $bla 1]
		} else {
		    set dnaseq [string range $data $from [expr {$to-1}]]
		    set qualseq [getqualfromto $dnaseq $quals $from $to]
                }
	
		if { $opts(-r) > 0 } {
		    if { [expr {$seqid % $opts(-r)}] == 0} {
			reverse_seq_qual dnaseq qualseq
		    }
		}

		set descline "1 [string length $dnaseq]"
		dumpfasta ${name}_$seqid $descline $dnaseq $fout $qualseq $qout
		
		# the very first and very last fragment will be put twice
		#  into the data set so that every part of the contig
		#  is covered at least twice
		if {$seqid == 1 || $lastflag >0 } {
		    # if we a re working with reversing sequences,
		    #  force a reverse of the sequence
		    if { $opts(-r) > 0 } {
			reverse_seq_qual dnaseq qualseq
		    }
		    dumpfasta ${name}_${seqid}_a $descline $dnaseq $fout $qualseq $qout
		}
	    }
        }
    }
}

proc fasta2frag::loadNextQuals {qin} {
    variable tempqualname

    set quals {}
    while {[gets $qin line] != -1} {
	if {[string compare [string index $line 0] ">"] == 0} {
	    set tempqualname "[string map {> ""} [lindex $line 0]]"
	    return $quals
	} else {
	    set quals [concat $quals $line]
	}
    }

    append quals $line

    return $quals
}

proc fasta2frag::processit {} {
    variable opts
    variable tempqualname

    set fin [open $opts(-infile) r]
    set fout [open $opts(-outfile) w]
    
    set havequal 0
    if {[file exist ${opts(-infile)}.qual]} {
	puts "have qual"
	set havequal 1
	set qin [open ${opts(-infile)}.qual r]
    } else {
	puts "no ${opts(-infile)}.qual"
	#set qout [open /dev/null w]
    }
    set qout [open ${opts(-outfile)}.qual w]

    set mpout ""
    set mbout ""
    if { $opts(-p) > 0 } {
	set mpout [open $opts(-outfile).pairs w]
	set mbout [open $opts(-outfile).bambus w]
	puts $mbout "library\tsimplepair\t[expr {$opts(-insert_size)-3*$opts(-insert_stdev)}]\t[expr {$opts(-insert_size)+3*$opts(-insert_stdev)}]"
    }

    gets $fin data
    set seqname "[string map {> ""} [lindex $data 0]]"

    if {$havequal} {
	gets $qin data
	set tempqualname "[string map {> ""} [lindex $data 0]]"
    }

    set actseq ""
    set cid 1

    while {[gets $fin line] != -1} {
	if {[string compare [string index $line 0] ">"] == 0} {
	    set quals {}
	    if {$havequal} {
		if {$tempqualname != $seqname} {
		    error "$seqname in fasta file is not equal to $tempqualname in quality file."
		}
		set quals [loadNextQuals $qin]
	    }
	    
	    fragASeq $actseq $seqname $fout $quals $qout $mpout $mbout
	    incr cid
	    set actseq ""
	    set seqname "[string map {> ""} [lindex $line 0]]"
	} else {
	    append actseq $line
	}
    }
    if {[string length $actseq] >0} {
	set quals {}
	if {$havequal} {
	    if {$tempqualname != $seqname} {
		error "$seqname in fasta file is not equal to $tempqualname in quality file."
	    }
	    set quals [loadNextQuals $qin]
	}
	fragASeq $actseq $seqname $fout $quals $qout $mpout $mbout
    }

    close $fin

    if { $opts(-p) > 0 } {
	close $mbout
	close $mpout
    }
}

proc usage {prgname} {
    puts stderr "$prgname: Splits a single fasta sequence into several
overlapping fragments.\n"
    puts stderr "Usage: $prgname ?options? infile outfile
\t-l   int      Length of fragments (default=3000)
\t-i   int      Increment of fragment start site (default=2500)
\t-p   int      Paired end (default=0 is off, 1 is on)
\t-s   int      Strobe sequencing (default=0 is off, 1 is on)
\t-q   int      Default quality when no quality data present (default=30)
\t-r   int      Reverse every nth fragment (not applicable with 
\t               paired-end where reverse one of the reads is reversed
\t               anyway (default=2)
\t-c   int      Circularise fragments so that they form a ring
                 (default=0 is is off, 1 would be on)

\t-qualdivisor  int      Divide quality values by this (default=1)
\t-minqual      int      But give it at least this qual (default=0)
\t-insert_size  int      paired-end: insert size (default=3000)
\t-insert_stdev int      paired-end: standard dev (default=900)
\t                        this is not working at the moment
\t-pairednaming string   naming scheme for paired-end: 
\t                        sanger, 454 (default) or solexa
\t-minmut       int      min. number of mutations/seq. errors (def=0)
\t-maxmut       int      max. number of mutations/seq. errors (def=0)
\t-strobeon     int      number of bases read during strobe on
\t-strobeoff    int      number of bases during strobe off
"
  exit
}

set fasta2frag::opts(-l) 3000
set fasta2frag::opts(-i) 2500
set fasta2frag::opts(-r) 2
set fasta2frag::opts(-p) 0
set fasta2frag::opts(-s) 0
set fasta2frag::opts(-q) 30
set fasta2frag::opts(-c) 0
#set fasta2frag::opts(-minmut) 1
#set fasta2frag::opts(-maxmut) 2
set fasta2frag::opts(-minmut) 0
set fasta2frag::opts(-maxmut) 0
set fasta2frag::opts(-pairednaming) 454
set fasta2frag::opts(-qualdivisor) 1
set fasta2frag::opts(-minqual) 0
set fasta2frag::opts(-insert_size) 3000
set fasta2frag::opts(-insert_stdev) 900
set fasta2frag::opts(-infile) ""
set fasta2frag::opts(-outfile) ""
set fasta2frag::opts(-strobeon) 100
set fasta2frag::opts(-strobeoff) 100



foreach {key val} $argv {
  if {![info exists fasta2frag::opts($key)]} {
      if {[string compare [string index $key 0] "-"] == 0} { 
	  puts stderr "Bad key $key\n"
	  usage $argv0
      }
      set val $key
      set key -infile
  }
  set fasta2frag::opts($key) $val
}

if {[string length $fasta2frag::opts(-infile)] ==0} { 
    puts "Missing '-infile' as keyword"
    usage $argv0 ;
}
if {[string length $fasta2frag::opts(-outfile)] ==0} { 
    puts "Missing '-outfile' as keyword"
    usage $argv0 ; 
}

fasta2frag::processit
