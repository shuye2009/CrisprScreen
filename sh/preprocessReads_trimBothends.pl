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

my $inFile = $ARGV[0];

# Global variables
# my $U6_str = substr("CGGTGTTTCGTCCTTTCCACAAGAT",0,$achorLen);    # Only use for 20 bases for now...should be enough
#my $TRACR_str = substr("GTTTTAGAGCTAGAAATAGCAAGTTAAAATA",0,$achorLen);
#my $pLCKO_TRACR_str = substr("CGGTACCTCGTCC",0,$achorLen);   # New stem sequence as per Oct 10, 2017
## read length is 31 in total, including VECTOR (stagger) length 0 to 5, guide length 20, anchor length 11 to 6
my $anchor_str = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATA"; # normally,  $achorLen = 11. For readLen = 31, guideLen = 20. This str will be adjusted based on read and guide
my $VECTOR_str = "CACCG"; # full length of stagger sequence
my $anchorLen = length($anchor_str);
my $mismatches = 0.15;  ## changed to fraction of anchor length by Shuye, allow 1 in 7 mismatch



# Routine for Hamming distance
sub hd {
   return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
 

sub match {

   if (hd($_[0],$_[1])/length($_[1]) <= $mismatches) {
      return 1;
   }

   return 0;
}
 
## handle split patterns, the stagger must match exactly, and the anchor can have mismatches
sub match2 {

   if ((hd($_[0],$_[1]) == 0 && (hd($_[2], $_[3])/length($_[3])) <= $mismatches)) {
      return 1;
   }

   return 0;
}


### or match partial anchor (length 6 to 10) exactly regardless of stagger
sub matchPartialAnchor {

   if (hd($_[0],$_[1]) == 0) {
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

# deprecated
sub searchForMatch {
   my $read = $_[0];
   my $anchor = $_[1];

   # The anchors are actually longer than what I'm using, so the furthest
   # 3' that the starting point can be would be the length of the anchor plus
   # a couple extra characters to account for what we're not using.  This solution
   # makes the function a little more generalized if the read length changes.
   my $len_index = 17; # original is 19, modified by Shuye to reduce the number of unmatched reads
   # $len_index is intrduced to allow guide length to be less than 20, 
   # min processed read length is 18 if len_index is set at 17, 
   # in the bowtie alignment step, max mismatch is 2, 
   # since the reference length is 20, reads less than 18 will treated as unliagned
   #
   for (my $i=(length($read)-length($anchor));$i>$len_index;$i--) {
      if (match( substr($read,$i,length($anchor)),$anchor)) {
         return $i;
      }
   }

   return -1;
}

# deprecated, subsumed by searchAtBothEnds
sub searchForMatchAt5prime {
   my $read = $_[0];
   my $vector = $_[1];

   # search for vector sequence at the 5 primr end. Handle the case where the start of vector is not the first nt
   my $len_index = 2; # only allow vector to start at first or second nt, not beyond second nt
   # $len_index is intrduced to allow guide length to be less than 20
   for (my $i=0;$i<$len_index;$i++) {
   	if (match( substr($read,$i,length($vector)), $vector)) {
   		return $i;
   	}
   }
 
   return -1;
}

# deprecated, not flexible
sub searchAtBothEnds_obsolete {
   my $read = $_[0];
   my @patterns = ("^CG\\w{20}GTTT\$", "^CCG\\w{20}GTT\$", "^ACCG\\w{20}GT\$");

   #search for partial vector sequence at the 5 primr end and partial anchor sequence at the 3 prime end.
   
   
   for (my $i=0;$i<scalar(@patterns);$i++) {
        if ($read =~ m/$patterns[$i]/) {
                   return $i+2;
        }
   }
   return -1;

}

sub searchAtBothEnds {
   my $read = $_[0];
   my $vector = $_[1];
   my $anchor = $_[2];

   #my @patterns = ("^CG\\w{20}GTTT\$", "^CCG\\w{20}GTT\$", "^ACCG\\w{20}GT\$");

   #search for partial or full vector sequence at the 5 primr end and partial anchor sequence at the 3 prime end.
   #start with sub-vector length 1 and full anchor, end with full vector, and sub-anchor length 6, 
   #allow guide length to go from 18 (two deletion) to 21 (one insertion)
   #allow vector starts at 1 only
   my $min_anchor = 3;
   my $best_diff = 100;  ## best_diff should be 0, as the real guide length is 20
   my $best = [-1, -1];
   for (my $len_subv=length($vector);$len_subv>=0;$len_subv--) {
   	for (my $len_guide=21;$len_guide>=18;$len_guide--) {
		#for (my $offset=0; $offset<=0; $offset++){ #allow zero random nt at the start of read, tried with one, but did not work well
  			 
	     	my $len_suba = length($read)-$len_subv-$len_guide;
	
		last if($len_suba < $min_anchor);  ## skip the read if remaining anchor length is less than $min_anchor

		my $suba = $anchor;
		if($len_suba < length($anchor)){
			$suba = substr($anchor, 0, $len_suba); # len_suba < anchorlen, set suba to substring of anchor
		}else{
			$len_suba = length($anchor); # if len_suba >= anchorlen, set it to anchorLen 
		}

		## if stagger exists, match stagger exactly and match anchor with mismatch, or match achor exactly
		if($len_subv > 0){
			my $subv = substr($vector, length($vector)-$len_subv, $len_subv);
        		if(match2(substr($read, 0, $len_subv), $subv, substr($read, length($read)-$len_suba, $len_suba), $suba) or matchPartialAnchor(substr($read, length($read)-$len_suba, $len_suba), $suba)){
				if(abs($len_guide - 20) < $best_diff){
					$best_diff = abs($len_guide - 20);
					$best = [$len_subv, $len_guide];
				}
			}
		## if stagger does not exist, match anchor with mismatch, as the anchor will be long enough
        	}elsif($len_subv == 0){
			if(match(substr($read, length($read)-$len_suba, $len_suba), $suba)){
				if(abs($len_guide - 20) < $best_diff){
                        		$best_diff = abs($len_guide - 20);
                                	$best = [$len_subv, $len_guide];
				}
                        }
		}
	
	}
   }

   if(@$best[1] != -1){
   	return $best;
   }else{
	return -1;
   }

}


my $anchor = $anchor_str;
my $vector = $VECTOR_str;
my %matchPos_stat = ();
my %matchLen_stat = ();
# Open files to start processing
open (IN,"$inFile") or die "Could not open input $inFile";
my $anchor_match = 0;
my $vector_match = 0;
my $both_match = 0;
my $unmatch = 0;
while (my $line1 = <IN>) {
#for (my $i=1;$i<15;$i++) {
#my $line1 = <IN>;
      chomp($line1);
      chomp(my $line2 = <IN>);
      chomp(my $line3 = <IN>);
      chomp(my $line4 = <IN>);


      my $guide = "";
      my $qual = "";

      # If 1 or 2 random nucleotides appear in front of the full vector, trim them
      # before search for guide

      for(my $offset=1; $offset<=3; $offset++){
		my $leader = substr($line2, $offset, length($vector));
      		if(hd($leader, $vector) == 0){
			$line2 = substr($line2, $offset);
			$line4 = substr($line4, $offset);
		}
	}

      # This solution is likely less efficient than using "match" to find the exact
      # string matchs, but this allows us to handle mismatches.  Keep searching for
      # more optimized solution, but once running on the server it doesn't matter
      # if it takes a while...
      my $matchPos = searchAtBothEnds($line2, $vector, $anchor);     
      if ($matchPos != -1) {
		my @pos = @$matchPos;
		my $len = $pos[1];
                my $start = $pos[0];
                $guide = substr($line2,$start,$len);
                $qual = substr($line4,$start,$len);
	
		if($start == 5){
			$vector_match++;
		}elsif($start == 0){
			$anchor_match++;
		}else{
			$both_match++;
		}

		if(exists $matchPos_stat{$start}){
                	my $count = $matchPos_stat{$start};
                	$count++;
                	$matchPos_stat{$start} = $count;
        	}else{
                	$matchPos_stat{$start} = 1;
        	}

		if(exists $matchLen_stat{$len}){
                        my $count = $matchLen_stat{$len};
                        $count++;
                        $matchLen_stat{$len} = $count;
                }else{
                        $matchLen_stat{$len} = 1;
                }

      } else {# Deal with cases where neither anchor is found (ie: unmatched)
       		$unmatch++;
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

print STDERR "position\tcount\n";
foreach my $pos (sort {$a <=> $b} keys %matchPos_stat){
	my $count = $matchPos_stat{$pos};
	print STDERR "$pos\t$count\n";
}
my $match_rate = 100*($anchor_match+$vector_match+$both_match)/($anchor_match+$vector_match+$both_match+$unmatch);
print STDERR "\nfull anchor $anchor_match; vector $vector_match; both $both_match\nunmatched $unmatch\nmatch rate $match_rate\nfinish\n";

print STDERR "guide\tcount\n";
foreach my $len (sort {$a <=> $b} keys %matchLen_stat){
        my $count = $matchLen_stat{$len};
        print STDERR "$len\t$count\n";
}

## END
