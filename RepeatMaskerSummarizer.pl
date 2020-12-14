#!/usr/bin/perl
use strict;
use warnings;

##Counts the number of ALL elements and number of nucleotides by both repeat type and class
##use:
# RepeatMaskerSummarizer.pl /path/to/repeatmasker.fa.out
# writes summary to /path/to/repeatmasker.fa.out_summary.txt

my $rmOutFile = $ARGV[0];
my $summaryFile = $rmOutFile ."_summary.txt";

my %repeatCount; 
my %repeatNtCount;
my %classCount; 
my %classNtCount;
my %familyCount; 
my %familyNtCount;
my $total = 0;
my $totalnt = 0;
		
open(RMOUT, "<", $rmOutFile);
open(SUMMARY, ">", $summaryFile);

while (my $line = <RMOUT>) {
	if($line =~ /^\s*SW   perc/){
		next;
	}
	if($line =~ /^\s*score/){
		next;
	}

	if($line =~ /^\s*$/){
		next;
	}

	$total +=1;
	chomp($line);
	$line =~ s/^ +//g;
	$line =~ s/ +/ /g;
 	my ($score, $percdiv, $percdel, $percins, $query, $qstart, $qend, $qleft, $dir, $repeat, $class, $rstart, $rend, $rleft, $id) = split(" ", $line);
	my $family = $class;			
	
	if($class =~ /\//){
		my @parts = split("\/", $class);
		$class = $parts[0];
		$family = $parts[1];
	}
				
	$repeatCount{$repeat} += 1;
	$classCount{$class} += 1;
	$familyCount{$family} += 1;
				
	my $qLen = $qend-$qstart + 1;
	$repeatNtCount{$repeat} += $qLen;
	$classNtCount{$class} += $qLen;
	$familyNtCount{$family} += $qLen;

	$totalnt += $qLen;
	
	#print SUMMARY "$class\t$repeat\t$repeatCount{$repeat}\t$repeatNtCount{$repeat}\t$classCount{$class}\t$classNtCount{$class}\n";
}
			
my @repeatkeys = keys(%repeatCount);
my @classkeys = keys(%classCount);

print SUMMARY "name\tbp\tcount\tbpPercent\tcountPercent\n";
print SUMMARY "\n";
print SUMMARY "............Class............\n";
my @printResults;
my $resPointer;
$resPointer = &printSummary(\%classNtCount, \%classCount, $totalnt, $total);
@printResults = @{$resPointer};
foreach my $line (@printResults) {
	print SUMMARY $line; 
}

print SUMMARY "............Family............\n";
$resPointer = &printSummary(\%familyNtCount, \%familyCount, $totalnt, $total);
@printResults = @{$resPointer};
foreach my $line (@printResults) {
	print SUMMARY $line; 
}

print SUMMARY "............Type............\n";
$resPointer = &printSummary(\%repeatNtCount, \%repeatCount, $totalnt, $total);
@printResults = @{$resPointer};
foreach my $line (@printResults) {
	print SUMMARY $line; 
}

print SUMMARY "............Total............\n";			
print SUMMARY "Total\t$total\t";
printf SUMMARY ("%.0f",$totalnt);
print SUMMARY "\n";
	

sub printSummary{
	my $bpPointer = shift;
	my $countPointer = shift;
	my $totalBP = shift;
	my $totalCount = shift;

	my %repBP = %{$bpPointer};
	my %repCount = %{$countPointer};
	
	my @results;
	my $dataLine;
	foreach my $key (sort keys %repBP) {
		my $percentTotalBP = $repBP{$key} / $totalBP;
		my $percentTotalCount = $repCount{$key} / $totalCount;
	
		$dataLine = "$key\t$repBP{$key}\t$repCount{$key}\t$percentTotalBP\t$percentTotalCount\n";

		push(@results, $dataLine);
	}
	return \@results;
}


#Once completed, go into excel.			
#Open the "summarize.txt" file in excel.
		 	


