#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2010 and later by Bastien Chevreux
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

namespace eval faqsel {
  variable opts

}


proc faqsel::conditionalFASTQCopy {fin fout} {
    variable names

    gets $fin line
    set seqname "[string range [lindex $line 0] 1 end]"
    set mustout 0
    if {[info exists names($seqname)]} {
	set seqname "[string range [lindex $line 0] 1 end]"
	set mustout 1
    }

    set seqdata ""
    set qualdata ""

    set inseq 1

    while {[gets $fin line] != -1} {
	#puts "l: $line"
	#puts "sn: $seqname"
	#puts "sd: $seqdata"
	#puts "qd: $qualdata"
	#puts "fl: $mustout $inseq"
	if {$inseq} {
	    if {[string compare [string index $line 0] "+"] == 0} {
		set inseq 0
		set sname "[string range [lindex $line 0] 1 end]"
		if {[string length sname] && $sname != $seqname} {
		    error "Last sequence name: $seqname\nNow reading line: $line\nThe names don't match?!"
		}
	    } elseif {[string compare [string index $line 0] "@"] == 0} {
		# sequence without qual
		if {$mustout} {
		    puts $fout "@$seqname"
		    puts $fout $seqdata
		}
		set seqdata ""
		set seqname "[string range [lindex $line 0] 1 end]"
		set mustout 0
		if {[info exists names($seqname)]} {
		    set mustout 1
		}
	    } else {
		append seqdata $line
	    }
	} else {
	    if {[string compare [string index $line 0] "@"] == 0} {
		if {[string length $seqdata] == [string length $qualdata]} {
		    if {$mustout} {
			puts $fout "@$seqname"
			puts $fout $seqdata
			puts $fout "+"
			puts $fout $qualdata
		    }
		    set seqdata ""
		    set qualdata ""
		    set seqname "[string range [lindex $line 0] 1 end]"
		    set inseq 1
		    set mustout 0
		    if {[info exists names($seqname)]} {
			set mustout 1
		    }
		} else {
		    append qualdata $line
		}
	    } else {
		append qualdata $line
	    }
	}
    }

    if {$mustout && [string length $seqdata]} {
	puts $fout "@$seqname"
	puts $fout $seqdata
	if {[string length $qualdata]} {
	    puts $fout "+"
	    puts $fout $qualdata
	}
    }
}

proc faqsel::processit {} {
    variable opts
    variable names

    puts "Reading names"
    set fin [open $opts(-name) r]
    while {[gets $fin line] != -1} {
#	set names([string trim $line]) 1
	set names([lindex $line 0]) 1
    }
    close $fin

    puts "Copying sequence data"
    set fin [open $opts(-infile) r]
    set fout [open $opts(-outfile) w]
    conditionalFASTQCopy $fin $fout
    close $fout
    close $fin
}

proc usage {prgname} {
    puts stderr "fastqselect: Select fastq sequences in a file according to names given 
             in a name file.\n"
    puts stderr "Usage: fastqselect ?options? 
\t-infile  name   filename containing all fasta sequences
\t-name    name   filename containing all names of sequences to select
\t-outfile name   filename where to write selcted sequences to
"
  exit
}

set faqsel::opts(-infile) ""
set faqsel::opts(-outfile) ""
set faqsel::opts(-name) ""

foreach {key val} $argv {
  if {![info exists faqsel::opts($key)]} {
      if {[string compare [string index $key 0] "-"] == 0} { 
	  puts stderr "Bad key $key\n"
	  usage $argv0
      }
      set val $key
      set key -infile
  }
  set faqsel::opts($key) $val
}

if {[string length $faqsel::opts(-infile)] ==0} { usage $argv0 ; }
if {[string length $faqsel::opts(-outfile)] ==0} { usage $argv0 ; }
if {[string length $faqsel::opts(-name)] ==0} { usage $argv0 ; }

faqsel::processit
