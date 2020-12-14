#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
    
my $faFile = $ARGV[0];
my $bedFile = $faFile .".bed";

my $in = Bio::SeqIO->new(-file   => "$faFile", -format => "fasta", );

open(BED, ">", $bedFile);

while (my $seq = $in->next_seq) {
    print BED $seq->id ."\t0\t". length($seq->seq) ."\n";
}
