#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2006 and later by Bastien Chevreux
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

namespace eval fassel {
  variable opts

}


proc fassel::conditionalFASTACopy {fin fout} {
    variable names

    gets $fin data
    set seqname "[string map {> ""} [lindex $data 0]]"
    set mustout 0
    if {[info exists names($seqname)]} {
	puts $fout ">$seqname"
	set mustout 1
    }

    while {[gets $fin line] != -1} {
	if {[string compare [string index $line 0] ">"] == 0} {

	    set seqname "[string map {> ""} [lindex $line 0]]"
	    if {[info exists names($seqname)]} {
		puts $fout ">$seqname"
		set mustout 1
	    } else {
		set mustout 0
	    }
	} elseif {$mustout >0} {
	    puts $fout $line
	}
    }
}

proc fassel::processit {} {
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
    conditionalFASTACopy $fin $fout
    close $fout
    close $fin

    if {[file exists ${opts(-infile)}.qual]} {
	puts "Copying quality data"
	set fin [open ${opts(-infile)}.qual r]
	set fout [open ${opts(-outfile)}.qual w]
	conditionalFASTACopy $fin $fout
	close $fout
	close $fin
    }

}

proc usage {prgname} {
    puts stderr "fastaselect: Select fasta sequences in a file according to names given 
             in a name file.\n
If fasta quality file is present (same basename, but with .qual appended),
then also selects sequences from there.\n"
    puts stderr "Usage: fastaselect ?options? 
\t-infile  name   filename containing all fasta sequences
\t-name    name   filename containing all names of sequences to select
\t-outfile name   filename where to write selcted sequences to
"
  exit
}

set fassel::opts(-infile) ""
set fassel::opts(-outfile) ""
set fassel::opts(-name) ""

foreach {key val} $argv {
  if {![info exists fassel::opts($key)]} {
      if {[string compare [string index $key 0] "-"] == 0} { 
	  puts stderr "Bad key $key\n"
	  usage $argv0
      }
      set val $key
      set key -infile
  }
  set fassel::opts($key) $val
}

if {[string length $fassel::opts(-infile)] ==0} { usage $argv0 ; }
if {[string length $fassel::opts(-outfile)] ==0} { usage $argv0 ; }
if {[string length $fassel::opts(-name)] ==0} { usage $argv0 ; }

fassel::processit
