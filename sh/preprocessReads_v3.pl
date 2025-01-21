#!/usr/bin/perl -w
#
# Need to write code to pre-process the new CRISPR v2.1 reads to
# allow for six possible read structures.  Either the tracr or the 
# reverse complement of the U6 promoter will be at the 3'-end of the 
# read, but the 5'-end can be padded by 1 to 8 bases.
#
# First version:  Strict string matching on N bases of tracr or U6 as anchor
# Second version:  Allow M mismatches in the anchor
#
# Date:  Sept 22, 2016
# Author:  Kevin R Brown
#
# INPUT:  FASTQ filename to parse
# OUTPUT:  Filtered FASTQ, Unmatched reads FASTQ
#
# ########################################

use strict;

## PUT IN SOME ERROR HANDLING TO CHECK FOR PROPER CALL TO FUNCTION
## CHECK NUMBER OF INPUT ARGUMENTS

my $achorLen = $ARGV[0];
my $inFile = $ARGV[1];

# Global variables
my $U6_str = substr("CGGTGTTTCGTCCTTTCCACAAGAT",0,$achorLen);    # Only use for 20 bases for now...should be enough
my $TRACR_str = substr("GTTTTAGAGCTAGAAATAGCAAGTTAAAATA",0,$achorLen);
my $pLCKO_TRACR_str = substr("CGGTACCTCGTCC",0,$achorLen);   # New stem sequence as per Oct 10, 2017
my $mismatches = 0.2;  ## changed to fraction of anchor length by Shuye, as the anchor length could be as short as 5 



# Routine for Hamming distance
sub hd {
   return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
 

sub match {

   if (hd($_[0],$_[1])/$achorLen <= $mismatches) {
      return 1;
   }

   return 0;
}
 

sub reverseComplement {
   my $seq = $_[0];
   chomp($seq);
   $seq =~ tr/ACGTacgt/TGCAtgca/;
   $seq=reverse($seq);
   return $seq;
}

sub searchForMatch {
   my $read = $_[0];
   my $anchor = $_[1];

   # The anchors are actually longer than what I'm using, so the furthest
   # 3' that the starting point can be would be the length of the anchor plus
   # a couple extra characters to account for what we're not using.  This solution
   # makes the function a little more generalized if the read length changes.
   my $len_index = 18; # original is 19, modified by Shuye to reduce the number of unmatched reads
   # $len_index is intrduced to allow guide length to be less than 20
   for (my $i=(length($read)-length($anchor));$i>$len_index;$i--) {
      if (match( substr($read,$i,length($anchor)),$anchor)) {
         return $i;
      }
   }

   return -1;
}

my $anchor = "";

# Open files to start processing
open (IN,"$inFile") or die "Could not open input $inFile";

while (my $line1 = <IN>) {
#for (my $i=1;$i<15;$i++) {
#my $line1 = <IN>;
      chomp($line1);
      chomp(my $line2 = <IN>);
      chomp(my $line3 = <IN>);
      chomp(my $line4 = <IN>);

      # Initialize the anchor
      if (length($anchor) == 0) {
         if ($line2 =~ m/$U6_str/g) {
            $anchor = $U6_str;
         } elsif ($line2 =~ m/$TRACR_str/g) {
            $anchor = $TRACR_str;
         } elsif ($line2 =~ m/$pLCKO_TRACR_str/g) {
            $anchor = $pLCKO_TRACR_str;
         }
         
         # Skip this block if not found
         if (length($anchor) == 0) {
            next;
         }else{
	    print "anchor $anchor\n";
	 }
      }

      my $guide = "";
      my $qual = "";

      # This solution is likely less efficient than using "match" to find the exact
      # string matchs, but this allows us to handle mismatches.  Keep searching for
      # more optimized solution, but once running on the server it doesn't matter
      # if it takes a while...
      if ((my $matchPos = searchForMatch($line2,$anchor)) > 0) { 
         my $len = 20;  ## this variable is added to handle the cases where the guide length is less than 20
	 my $start = $matchPos-$len;
	 if($start < 0){
		$len += $start;
		$start = 0;
		
	 }
		

         if ($anchor eq $U6_str) {
            $guide = reverseComplement(substr($line2,$start,$len));
            $qual = reverse(substr($line4,$start,$len));
         } else {
            $guide = substr($line2,$start,$len); 
            $qual = substr($line4,$start,$len); 
         }
      } else {      # Deal with cases where neither anchor is found (ie: unmatched)
         print STDERR $line1."\n";
         print STDERR $line2."\n";
         print STDERR $line3."\n";
         print STDERR $line4."\n";
         next;
      }

      print STDOUT $line1."\n";
      print STDOUT $guide."\n";
      print STDOUT $line3."\n";
      print STDOUT $qual."\n";

} # End of while


# Close all the files
close(IN);
 
