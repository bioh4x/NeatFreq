This is an artificial set for ecoli (first 5040 bases) just to show a few
things.

This simulates mapping of a genome (strain named "test") against the first
5040 bases of ecoli K12 MG1655 (which were circularised as if it were a full
genome). Simply execute the runme.sh script.

MIRA automatically tags bases of interest. If you convert the CAF result file
to a gap4 database (caf2gap -project test -ace ecoli_out.caf), you will see

1) a deletion of 6 bases around position 255. Indels of this size start to
   get large for finding them with Solexa 36mers, MIRA had trouble aligning
   some of the reads (you might want to realign a bit by hand). But the
   location is plentily advertised by MIRA with a SROc tag and WRMc tags 
2) a clean base change at position ~362 (tagged SROc in the consensus and SROr
   in the reads)
3) a good one base indel at position ~482 (allso tagged SROc)
4) a 3+1 base insertion at position ~617 which also gave MIRA a bit of
   trouble, but again the location is advertised by SROc and WRMc tags. When
   working on a cleanup in gap4, I'd realign this area quickly by hand, and
   delete the WRMc tags.



