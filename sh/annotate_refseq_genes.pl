#!/usr/bin/perl -w
#
## this is a helper function for providing gene information for refseq gene ids
## to be used in process_CIMS_CITS_results.sh pipeline
#
## also produce a filtered bed file that only contain sites located in 3' UTR
my $TARGET = $ARGV[0]; ## name of peaks
my $input_file = $ARGV[1]; # an output from homer annotatePeak.pl: ${TARGET}_CITS_crosslinkSite_${motif1}_annotation.txt 
my $wd = $ARGV[2];
chdir $wd;

print "Processing $input_file\n";

my $HOME = $ENV{HOME};
my $annotation_file = $HOME."/homer/data/accession/human2gene.tsv";
my $gene_file = $HOME."/homer/data/accession/human.description";
my $output_file = $input_file;
$output_file =~ s/\.txt/_description\.txt/g;
my $filter_bed_file = $input_file;
my %refseq_description_map = ();
my %refseq_gene_map = ();

open(D, "<$annotation_file") or die "$annotation_file error!";
open(G, "<$gene_file") or die "$gene_file error!";

while(<D>){
	my $tmp = $_;
	chomp $tmp;
	$tmp =~ s/\n//;
	$tmp =~ s/\r//;

	my @tary = split(/\t/, $tmp);
	if(defined $tary[0]){
		my $refseq = $tary[0];
		$refseq_description_map{$refseq} = $tmp;
	}
}
close(D);

while(<G>){
        my $tmp = $_;
        chomp $tmp;
        $tmp =~ s/\n//;
        $tmp =~ s/\r//;

        my @tary = split(/\t/, $tmp);
        if(defined $tary[2]){
                my $refseq = $tary[2];
                $refseq_gene_map{$refseq} = $tary[9]."\t".$tary[8];
        }
}
close(G);

open(I, "<$input_file");
open(O, ">$output_file");

my $ih = <I>;
$ih =~s/\n//;
print O "Chr\tStart\tEnd\tFeature\tScore\tStrand\tID\tEntrez\tUnigene\tRefseq\tEnsembl\tNote\tGeneName\tGeneType\tDescription\n";
while(<I>){
        my $tmp = $_;
        chomp $tmp;
        $tmp =~ s/\n//;
        $tmp =~ s/\r//;

        my @tary = split(/\t/, $tmp);
        if(defined $tary[7]){
		my $refseq_info = $tary[7];
                #print "$refseq_info\n";

		# Homer extends the peak to a given size, find the original center of peak
		my $start = $tary[2] + int(($tary[3]-$tary[2])/2) - 1;
                my $end = $start + 1;
		my $bedline = "$tary[1]\t$start\t$end\t$refseq_info\t$tary[5]\t$tary[4]";

		if($refseq_info =~ /.*\(([N|X][M|R]_\d+).*\)/){
			my $accession = $1;
			#print "$accession\n";
			if(defined $refseq_description_map{$accession}){
				
				my $description = $refseq_description_map{$accession};
				my @des = split(/\t/, $description);
				my $refseq = $des[3];
				my $gene = $refseq_gene_map{$refseq};
        		        print O "$bedline\t$description\t$gene\n";
	
			}else{
				print O "$bedline\t\t\t\t\t\t\t\tOther\t\n"; # for transcripts that are obsolete (reasigned or removed)
			}
		}elsif($refseq_info =~ /^Intergenic$/){
			print O "$bedline\t\t\t\t\t\t\t\tIntergenic\t\n";
			#print "no matching found for $tmp\n";
		}

	}

}
close(I);
close(O);
print "finish adding extra gene information\n";
