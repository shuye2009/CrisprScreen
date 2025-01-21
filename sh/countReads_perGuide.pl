#!/usr/bin/perl -w
#
# organize read count data by gene
#
# Date:  Nov 15, 2019
# Author:  Shuye Pu
# Updated: Mar 06, 2021
#
# INPUT:  two column count data: count\tguide
# OUTPUT:  Nine columns count data: gene\tguide1\tcount1\tguide2\tcount2\tguide3\tcout3\tguide4\tcount4
#
# ########################################

use strict;


my $inFile = $ARGV[0];
my $outFile = $ARGV[1];
my $scaled = $ARGV[2]; ## take T for true F for false
my $total_readCount = 0;

# Global variables
my $HOME = $ENV{HOME};
my $guide_library = "$HOME/Zuyao/crispr/TKO/TKOv2.1-Human-Library.txt";
my %gene_guide_count = ();
my %guide_sequence = ();
open(G, "<$guide_library") or die "Could not open library $guide_library";
my $header = <G>;

while (my $line = <G>) {
      chomp($line);
      my @tary = split(/\t/, $line);
      my $gene = $tary[2];
      my $guide = $tary[1];
      my $sequence = $tary[0];
      $guide_sequence{$guide} = $sequence;
      $gene_guide_count{$gene}{$guide} = 0;
}

# Open files to start processing
open (IN,"<$inFile") or die "Could not open input $inFile";
open (OUT, ">$outFile");

while (my $line = <IN>) {

      chomp($line);
      my @tary = split(/\s+/, $line);
      my $count = $tary[0];
      my $guide = $tary[1];
      if(scalar(@tary) == 3){
      	$count = $tary[1];
      	$guide = $tary[2]; 
      }
      #print "$line\n";
      #print "$guide\n";
      my @parts = split(/_/, $guide);
      my $gene = $parts[1];
      my $sequence = $guide_sequence{$guide};
 
      print OUT "$sequence\t$gene\t$count\n";

      #print "$gene\n";
      $gene_guide_count{$gene}{$guide} = $count;
      $total_readCount += $count;

} # End of while
# Close all the files
close(IN);
 
if($scaled eq "T"){
	$total_readCount /= 10000000; ## scale to ten millions
}else{
	$total_readCount = 1;
}
print STDERR "Total read count scale factor: $total_readCount\n";
# produce output
print STDOUT "Guide1\tCount1\tGuide2\tCount2\tGuide3\tCount3\tGuide4\tCount4\n";
my @genes = keys %gene_guide_count;
foreach my $gene (sort(@genes)){
	my $ref = $gene_guide_count{$gene};
	my %guide_count = %$ref;
	my @guides = sort(keys(%guide_count)); ## to keep the order of guides consistent over files
	print STDOUT "$gene";
	
	for(my $i=0; $i<4; $i++){
		my $guide = "NA";
		my $count = "NA";
		if(defined $guides[$i]){
			$count = $gene_guide_count{$gene}{$guides[$i]}/$total_readCount;
			$count = sprintf("%.1f", $count);
			$guide = $guides[$i];
		}
		print STDOUT "\t$guide\t$count";
	}
	print STDOUT "\n";
		
}
