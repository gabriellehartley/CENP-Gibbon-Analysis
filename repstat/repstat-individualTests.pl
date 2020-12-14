#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#quick and easy way to parse all of the signfican outputs, you might need to tweak where described

#grep -v Bed1.counts comparison/cdfPlots/summaryStats---*.csv | sed "s/\"//g" | sed "s/,/\t/g" | awk '$3 > 10 || $4 > 10' > temp && \ #the $3 and $4 sets some minimums for the counts
#TESTS=`wc -l < temp |xargs echo -n` && cat temp | awk -v var="$TESTS" '{print $0 ,"\t", ( 0.01 / var)}' | awk '($7 + 0) < $9' | \ #must be above bonferroni corrected at 0.01 / var
#sort -k7,7g | grep "\-species\-" | cut -f2 | sort -u | tr "\n" "|" | sed 's/|/\\|/g' | sed 's/.\{2\}$/\n/' | \
#while read -r REPLACE; do grep $REPLACE fasta-RM.out > sig-comparison.out; done && rm temp #fasta-RM.out is the original rm.out file

my $versionNumber = "0.11";

my $rmOutFile;
my $bed1File;
my $bed2File;
my $sSize;
my $wSize;
my $wd;
my $cores;
my $minus21 = 0;
my $skipR;
my $runR = 1;
my $clean = 1;
my $debug = 0;
my $runITests = 0;
my $runClustTests = 0;
my $runClassTests = 0;
my $justR = 0;
my $runAll = 0;
my $version = 0;
my $mcLength = 100000;
my $multiLength = 1;
#c("percent_divergence", "percent_deleted", "percent_inserted", "length", "lenRegRepNorm", "lenRegDNANorm", "lenAllRepNorm", "lenAllDNANorm", "percent_identity", "bpStart", "bpStop")
my $runIMetrics = "T,T,T,T,T,T,T,T,T,T,T";

my $runCMetrics = "T,T,T,T,T,T,T,T,T";
my $runGTTests;
my $runLTTests;
my $run2STests;
my $help;

GetOptions ("rmOut=s"       => \$rmOutFile,
            "bed1=s"        => \$bed1File,
	    "bed2=s"        => \$bed2File,
	    "minus2_1"      => \$minus21,
	    "outDir=s"      => \$wd,
	    "cores=i"       => \$cores,
	    "noR"           => \$skipR,
	    "justR"         => \$justR,
	    "clean"         => \$clean,
	    "debug"         => \$debug,
	    "individual"    => \$runITests,
	    "clustering"    => \$runClustTests,
	    "classify"      => \$runClassTests,
	    "runAll"        => \$runAll,
	    "version"       => \$version,
            "help"          => \$help,
	    "runCMetrics"   => \$runCMetrics,
            "runIMetrics"   => \$runIMetrics,
	    "runGTTests"    => \$runGTTests,
	    "runLTTests"    => \$runLTTests,
	    "run2STests"    => \$run2STests,
            "mcLength"      => \$mcLength )
    or die("Error in command line arguments, use --help for more information\n");


if($version){
    print "repStat version: $versionNumber \n";
    exit;
}


if($help || !($rmOutFile && $bed1File && $bed2File && $wd && $cores)){
    print "\n";
    print "Usage:\n";
    print "repStat --cores 10 --rmOut seq.fa.out --bed1 ROI.bed --bed2 background.bed --outDir outputDirectory\n";
    print "\n";
    print "Options are:\n";
    print "--rmOut \'file\'           repeatmasker .out file\n"; #extra spaces are because the \' takes up only one char when printed
    print "--bed1 \'file\'            a bed file defines the region of interest\n";
    print "--bed2 \'file\'            a bed file defines a second region of interest / background region\n";
    #print "--minus2_1               if the repeats found in bed1 should be removed from bed2\n";
    print "--outDir \'directory\'     where the output should be saved\n";
    print "--cores \'num\'            number of cores for R to use, the second portion of the program maxes out at 10 cores\n";
    print "--noR                    disables the R portion of the program\n";
    print "--justR                  runs just the R portion of the program\n";
    #print "--runAll                 tells repStat to perform a full analysis of the repeats \(default\, overrides --individual\, --clustering\ and --classify\)\n";
    #print "--individual             tells repStat to perform the tests for individual repeats \(Wilcox Signed Rank, Mann Whitney U, CDFs\, violins and summary stats\)\n";
    #print "--mcLength               tells the coin library how many samples should be taken for the Monte Carlo estimation of the p-values, defaults to 100,000\n";
    print "--runIMetrics            takes a list of T or F to pick which tests of individual repeats, see details for specifics\n";
    print "--runGTTests             this will perform one sided tests for if the parameter tested in bed1 is greater than in bed2\n";
    print "--runLTTests             this will perform one sided tests for if the parameter tested in bed1 is less than in bed2\n";
    print "--run2STests             this will perform two sided tests for if the parameter tested in bed1 is simply different from in bed2, only needs to be set explicitly if another is set\n";
    #print "--clustering             tells repStat to performs the clustering analyses for the repeats \(MDS\/PCA\, heatmaps\, pvclust and pie charts\)\n";
    #print "--classify               tells repStat to performs the classification analyses for the repeats \(knn, kMeans\, MCLUST\)\n";
    #print "--tidy                   compresses all of the raw data\n";
    print "--clean                  deletes all of the files R needed to generate the figures\n";
    print "--version                prints the version number and then quits\n";
    print "--help \n";
    print "\n";
    print 'repStat depends on R having: "gplots", "parallel", "ggplot2", "Matching", "coin", "data.table", "ape"'. "\n";
    print "and the PATH variable to include bedtools\n";
    print "\n";
    print "Running Suggestions:\n";
    print "Use the p-values from the outDir\/individualRepeatTests\/combinedStats--all--\*\.csv files to ID the repeats that are actually significantly different.\n";
    print "Then use grep to remove the uninteresting repeats from the rm.out file \(look at the \-v to invert the match\) before running the clustering step.\n";
    print "This can be done using: \"grep \'pattern1\\\|pattern2\' sequenceFile.fa.out\", be certain to include a space to avoid subpattern matching\n";
    print "This will make the MDS plots and pie charts generated by R a lot clearer.\n";
    print "\n";
    print "It's recommended to run repStat at least twice on the same sample and then compare the p-values to see if \"coin\" has converged on an\n";
    print "accurate measurement for the p-values, which can be a problem with infrequently occuring repeats.\n";
    print "\n";
    print "Details:";
    print "As it stands, repStat uses Holms to adjust the p-valeus to avoid false positive errors that arise from doing a large number of tests. To\n";
    print "help increase the statistical power of repStat use the \"--runIMetrics option to reduce the number of statistical tests performed.\" \n";
    print "--runIMetrics takes an array of 11 \"T\"s or \"F\"s corresponding to which tests to perform,  so to tests all 9 metrics, the following would be\n";
    print "used \"--runIMetrics T,T,T,T,T,T,T,T,T,T,T\", which is the default. The following is a list of all of the metrics and the abbrivation used in corresponding order:\n";
    print "1. percent divergence, 2. percent deleted, 3. percent inserted, 4. length of repeat,\n";
    print "5. a repeat's length as a fraction of the interval's repetitive DNA, 6. a repeat's length as a fraction of the interval's length,\n";
    print "7. a repeat's length as fraction of the total repetitive DNA, 8. a repeat's length as a fraction of the total DNA\n";
    print "9. percent identity, 10. first bp of repeat for region, 11. last bp of repeat for region\n";
    print "\n";
    print "1. percent_divergence, 2. percent_deleted, 3. percent_inserted, 4. length, 5. lenRegRepNorm, 6. lenRegDNANorm, 7. lenAllRepNorm, 8. lenAllDNANorm, 9. percent_identity, 10. bpStart, 11. bpStop\n";
    print "Please note that 1-4,10,11 are directly from the RepeatMasker out file and describe the aligned portion of the repeat\n";
    print "for 5 and 7, RepeatMasker can call the same base as belonging to two different repeats, so this is calculated by the amount of DNA that is called as a repeat\n";
    print "for 6 and 8, this is done using the total length of the supplied bed files, not overlapping";
    print "\n";
    print "Strongly consider reducing the number of metrics tested and which categories of repeats are of interest \*before\* running to minimize p-hacking\n";
    print "\n";
    print "\n";
    exit;
}

my($filename1, $dirs1, $suffix1) = fileparse($bed1File);
my($filename2, $dirs2, $suffix2) = fileparse($bed2File);
$filename1 =~ s/\.bed$//;
$filename2 =~ s/\.bed$//;

if($filename1 eq $filename2){
    $filename1 = $filename1 ."-bed1";
    $filename2 = $filename2 ."-bed2";
}


if($skipR){
    $runR = 0;
}

if($runITests || $runClustTests || $runClassTests){
    #at least one option has been set
}else{
    $runClassTests = 1;
    $runClustTests = 1;
    $runITests = 1;
}


if(defined($runGTTests) || defined($runLTTests)){
    if(defined($run2STests)){
	$run2STests = 1;
    }else{
	$run2STests = 0;
    }
}else{
    $run2STests = 1;
    $runGTTests = 0;
    $runLTTests = 0;
}


if (-e $bed1File && -e $bed2File && -e $rmOutFile){
    print "all files exist, ready to go\n";
}else{
    print "one of the supplied files doesn't exist\n";
    exit;
}

#my $wd = "workingFull1";
$wd =~ s/\/$//; #remove the trailing / in case it's supplied
system("mkdir -p $wd");

#   SW   perc perc perc  query         position in query              matching          repeat              position in repeat
#score   div. del. ins.  sequence      begin    end          (left)   repeat            class/family      begin  end    (left)       ID

# 7361   11.0  1.8  0.3  Contig0              3     1146    (44619) C rnd-2_family-58   LINE/L1             (55)   5378   4218        1
#  455   20.5  7.7  4.6  Contig0           1150     1318    (44447) C rnd-5_family-1113 Unknown             (91)    450    277        2
#  851   18.7  6.3  3.8  Contig0           1369     1574    (44191) C rnd-4_family-707  Unknown              (1)    218      8        3
#  730   19.8  3.4  9.7  Contig0           1937     2214    (43551) C rnd-4_family-1381 Unknown             (55)    342      5        4

if(!$justR){
    open(RMOUT, "<", $rmOutFile);

    print "loading in rm out file\n";
    my @file = <RMOUT>;
    close(RMOUT);
    #create bed file for bedtools to intersect with
    #while doing that read in the information to be printed
    my %repData;
    my %id2c;
    my %id2f;
    my %id2rep;

    my %sList;
    my %cList;
    my %fList;

    my %repOnScaf;

    open(BOUT, ">", $wd ."/allRM.bed");
    print "processing rm out file\n";
    #repeatmasker sometimes splits repeats that have the same ID number into multiple segments... %&!@$
    #so now the identifier in the code is the "id_seg_scaffold_at_startPosition" to allow for multiple repeatmasker out files to be merged
    #also the repeats can overlap, so careful about that
    foreach my $line (@file){
	#the following if series handles multiple repeatmasker files joined
	if($line =~ /^   SW/){
	    next;
	}elsif($line =~ /^score/){
	    next;
	}elsif($line =~ /^$/){
	    next;
	}

	chomp($line);

	$line =~ s/^\s+//g;
	#print $line ."\n";
	$line =~ s/\s+/\t/g;
	my @data = split(/\t/, $line);
	#collect the data
	my $id = $data[14];
	if(!defined($id)){
	    print "Skipping odd line:\n". $line ."\n";
	    next;
	}

	$id = $id ."_seg_" . $data[4] ."_at_". $data[5];
	$data[9] =~ s/\//_/g;
	$id2rep{$id} = $data[9];

	$data[10] =~ s/\?//;#fail
	if($data[10] =~ /\//){
	    my ($c,$f) = split(/\//, $data[10]);
	    $id2c{$id} = $c;
	    $id2f{$id} = $f;

	    $cList{$c} = 1;
	    $fList{$f} = 1;
	}else{
	    #if the repeat doesn't have a family, count the class as it's family
	    $id2c{$id} = $data[10];
	    $cList{$data[10]} = 1;
	    $id2f{$id} = $data[10];
	    $fList{$data[10]} = 1;
	}

	$sList{$data[9]} = 1;

	$id2rep{$id} = $data[9];
	$data[6]++;		#because bed is 0-based half open
	my $length = $data[6] - $data[5];

	my $startInRep = $data[11];
	$startInRep =~ s/\(//;
	$startInRep =~ s/\)//;
	my $stopInRep = $data[12];
	$stopInRep =~ s/\(//;
        $stopInRep =~ s/\)//;

	#my results = contig div% del% ins% length start stop
	my $results = $data[4] ."\t". $data[1] ."\t". $data[2] ."\t". $data[3] ."\t". $length ."\t". $startInRep ."\t". $stopInRep;
	$repData{$id} = $results;

	$repOnScaf{$data[4]} = 1;

	print BOUT "$data[4]\t$data[5]\t$data[6]\t$id\n";
    }
    @file = ();
    close BOUT;

    my $sortRules = "-T ". $wd ."/ -k1,1 -k2,2n";
    
    my $repSFile = $wd ."/allRepSorted.bed";
    my $cmd = "sort --parallel=". $cores ." $sortRules ". $wd ."/allRM.bed \> $repSFile";
    #print "$cmd\n";
    system($cmd);

    print "categorizing repeats based on supplied bed files and collecting statistics\n";

    #sorting the files
    my $b1s = $wd ."/bed1Sorted.bed";
    $cmd = "sort --parallel=". $cores ." $sortRules ". $bed1File ." \| cut -f1\,2\,3 \> $b1s";
    #print "$cmd\n";
    system($cmd);

    my $b2s = $wd ."/bed2Sorted.bed";
    $cmd = "sort --parallel=". $cores ." $sortRules ". $bed2File ." \| cut -f1\,2\,3 \> $b2s";
    #print "$cmd\n";
    system($cmd);

    my $fakeGenomeFile = $wd ."/constructed.genome";
    $cmd = "cat $repSFile $b1s $b2s | sort -T ". $wd ."/ -k1,1 -k3,3nr | awk \'\!_\[\$1\]\+\+\' | cut -f1,3 | awk \'\$2\+\=1\' | sed \'s\/ \/\\t\/\' > $fakeGenomeFile";
    #print "$cmd\n";
    system($cmd);

    #get the repeats for each region
    my $bed1RepFile = $wd ."/bed1Rep.bed";
    my $bed2RepFile = $wd ."/bed2Rep.bed";

    $cmd = "bedtools intersect -g $fakeGenomeFile -nonamecheck -wo -a ". $repSFile ." -b $b1s | sort $sortRules > $bed1RepFile\n";
    #print $cmd;
    system($cmd);

    $cmd = "bedtools intersect -g $fakeGenomeFile -nonamecheck -wo -a ". $repSFile ." -b $b2s | sort $sortRules > $bed2RepFile\n";
    #print $cmd;
    system($cmd);

    if($minus21){  #remove any repeat that was found in 1 from 2 by...
	#...first removing them from the general pool of repeats
	$cmd = "bedtools subtract -g $fakeGenomeFile -nonamecheck -A -a $repSFile -b $bed1RepFile > ". $wd ."/allRepeats-minus-bed1.bed\n";
	#print $cmd;
	system($cmd);
	#...then by intersecting that depleted repeat file
	$cmd = "bedtools intersect -g $fakeGenomeFile -nonamecheck -wo -a ". $wd ."/allRepeats-minus-bed1.bed -b $b2s \| sort $sortRules > $bed2RepFile\n";
	#print $cmd;
	system($cmd);

    } #doing it in this way lets the bed2 intervals remain intact, useful for later processing

    print "reading in the intervals for the bed files to get lengths of intervals\n";
    my $bed1DNALen;
    my $bed2DNALen;

    my %size1;
    open(BED1, $b1s);
    foreach my $line (<BED1>){
	chomp $line;
	my @data = split(/\t/, $line);
	my $name = $data[0] .":". $data[1] ."-". $data[2];
	$size1{$name} = $data[2] - $data[1];
	$bed1DNALen += $data[2] - $data[1];
    }
    close(BED1);

    my %size2;
    if($minus21){
	$bed2DNALen = 0;
	my $bed2Minus = $wd ."/bed2-1.bed";
	$cmd = "bedtools intersect -g $fakeGenomeFile -nonamecheck -wao -sorted -a $b1s -b $b2s > $bed2Minus";
	#print "$cmd\n";
	system($cmd);

	open(BED, $bed2Minus);

	foreach my $line (<BED>){
	    chomp $line;
	    my @data = split(/\t/, $line);
	    if($data[4] != "-1"){
		my $len = $data[5] - $data[4] - $data[6];
		$bed2DNALen += $len;
		my $interval = $data[3] .":". $data[4] ."-". $data[5];
		$size2{$interval} += $len;
	    }
	}
	close(BED);
    }else{
	open(BED2, $b2s);
	foreach my $line (<BED2>){
	    chomp $line;
	    my @data = split(/\t/, $line);
	    my $name = $data[0] .":". $data[1] ."-". $data[2];
	    $size2{$name} = $data[2] - $data[1];
	    $bed2DNALen += $data[2] - $data[1];
	}
	close(BED2);
    }

    my @bed1ID;
    my @bed2ID;

    print "processing intersect bed files\n";
    my %b1IntervalRepLen;
    my %b2IntervalRepLen;

    #get the amount of repeat dna per interval
    my $flattenedReps = $wd ."/repsFlattened.bed";
    $cmd = "bedtools merge -nonamecheck -i $repSFile > $flattenedReps";
    #print "$cmd\n";
    system($cmd);

    #the merging has to be done because sometimes repeatmasker says the same dna belongs to two different repeats
    my $bed1RepFlat = $wd ."/bed1RepFlat.bed";
    $cmd ="bedtools merge -nonamecheck -i $bed1RepFile | sort $sortRules | bedtools intersect -g $fakeGenomeFile -nonamecheck -sorted -wao -a $b1s -b stdin > $bed1RepFlat";
    #print $cmd ."\n";
    system($cmd);

    my $bed2RepFlat = $wd ."/bed2RepFlat.bed";
    $cmd ="bedtools merge -nonamecheck -i $bed2RepFile | sort $sortRules | bedtools intersect -g $fakeGenomeFile -nonamecheck -sorted -wao -a $b2s -b stdin > $bed2RepFlat";
    #print $cmd ."\n";
    system($cmd);

    my $bed1RepLen; 
    my $bed2RepLen; 

    open( BED1F, $bed1RepFlat);
    foreach my $line (<BED1F>){
	chomp($line);
	my @data = split(/\t/, $line);
	my $interval = $data[0] .":". $data[1] ."-". $data[2];
	$b1IntervalRepLen{$interval} += $data[6];
	#print "$b1IntervalRepLen{$interval} $interval $data[6]\n";
	$bed1RepLen += $data[6];
    }

    open( BED2F, $bed2RepFlat);
    foreach my $line (<BED2F>){
	chomp($line);
	my @data = split(/\t/, $line);
	my $interval = $data[0] .":". $data[1] ."-". $data[2];
	$b2IntervalRepLen{$interval} += $data[6];
	$bed2RepLen += $data[6];
    }

    #associate individual repeats to specific intervals (which can be overlapping)
    my $bed1RepWAO = $wd ."/bed1RepWAO.bed";
    $cmd ="bedtools intersect -g $fakeGenomeFile -nonamecheck -sorted -wao -a $b1s -b $bed1RepFile | cut -f1,2,3,4,5,6,7,12 | uniq > $bed1RepWAO";
    #print $cmd ."\n";
    system($cmd);

    my $bed2RepWAO = $wd ."/bed2RepWAO.bed";
    $cmd ="bedtools intersect -g $fakeGenomeFile -nonamecheck -sorted -wao -a $b2s -b $bed2RepFile | cut -f1,2,3,4,5,6,7,12 | uniq > $bed2RepWAO";
    #print $cmd ."\n";
    system($cmd);

    my %repSizes1;
    my %repCountsInterval1;
    my $bed1RepCount;

    open(BED1WAO, $bed1RepWAO);
    foreach my $line (<BED1WAO>){
	chomp($line);
	my @data = split(/\t/, $line);
	my $interval = $data[0] .":". $data[1] ."-". $data[2];
	my $rep = $data[6];
	if($data[4] != -1){
	    my $id = $rep ."_on_". $interval;
	    push(@bed1ID, $id);
	    $repSizes1{$id} = $data[7];
	    $repCountsInterval1{$interval}++;
	    $bed1RepCount++;
	}else{			  
	    next;
	}
    }
    close(BED1WAO);

    my %repSizes2;
    my %repCountsInterval2;
    my $bed2RepCount;

    open(BED2WAO, $bed2RepWAO);
    foreach my $line (<BED2WAO>){
	chomp($line);
	my @data = split(/\t/, $line);
	my $interval = $data[0] .":". $data[1] ."-". $data[2];
	my $rep = $data[6];
	if($data[4] != -1){
	    my $id = $rep ."_on_". $interval;
	    push(@bed2ID, $id);
	    $repSizes2{$id} = $data[7];
	    $repCountsInterval2{$interval}++;
	    $bed2RepCount++;
	}else{			  
	    next;
	}
    }
    close(BED2WAO);

    #write the data files

    print "writing results to files\n";
    #for the individual repeats

    my @categories = ("species", "class", "family");
    my @hashCon = (\%id2rep, \%id2c, \%id2f);

    my %avgDivRepSeq;
    my %avgDelRepSeq;
    my %avgInsRepSeq;
    my %avgIdtRepSeq;
    my %avgLenRepSeq;
    my %avgStopRepSeq;
    my %avgStartRepSeq;

    my %sumRepLenSeq;
    my %sumRepCountSeq;

    #these are useful if repeat changes bp quantity in sample
    my %sumRepLenNormRepTotalLen;
    my %sumRepLenNormDNATotalLen;
    my %sumRepLenNormRepIntervalLen;
    my %sumRepLenNormDNAIntervalLen;

    my %sumRepCountNormRepTotalLen;
    my %sumRepCountNormDNATotalLen;
    my %sumRepCountNormRepIntervalLen;
    my %sumRepCountNormDNAIntervalLen;

    my %sumRepCountNormRepTotalCount;
    my %sumRepCountNormDNATotalCount;
    my %sumRepCountNormRepIntervalCount;
    my %sumRepCountNormDNAIntervalCount;

    my %sumRepCountNormRepTotal;
    my %sumNormDNARepCount;

    my %b1Intervals;
    my %b2Intervals;

    #for pie charts
    my %pieCount;
    my %pieLength;

    #typeOfRep  contig/interval  sizeOfInterval(minusROI)
    my $b1StatsFile = $wd . "/bed1Stats.tsv";
    open(B1DATA, ">", $b1StatsFile);
    my $b2StatsFile = $wd . "/bed2Stats.tsv";
    open(B2DATA, ">", $b2StatsFile);

    my %sUsed;
    my %cUsed;
    my %fUsed;

    my $totalStats = $bed1DNALen ."\t". $bed2DNALen ."\t". $bed1RepLen ."\t". $bed2RepLen ."\t". $bed1RepCount ."\t". $bed2RepCount;

    for(my $i = 0; $i < @categories; $i++){
	my $c = $categories[$i];
	my %id2cat = %{$hashCon[$i]};
	print "on $c\n";
	#print "on bed1 regions\n";
	#for the family,but not every repeat has a family
	
	my $totalStats1 = $bed1DNALen ."\t". $bed2DNALen ."\t". $bed1RepLen ."\t". $bed2RepLen ."\t". $bed1RepCount ."\t". $bed2RepCount;
	foreach my $repID (@bed1ID){
	    my @data = split(/_on_/, $repID);
	    my $interval = $data[1];
	    my $id = $data[0];
	    my $contigBedRepLen = $b1IntervalRepLen{$interval};
	    
	    my $info;
	    my $whichBed = "bed1";
	    my @rData = split(/\t/, $repData{$id});
	    my $repDiv = $rData[1];
	    my $repDel = $rData[2];
	    my $repIns = $rData[3];
	    my $pIdentity = 100 - ($repDiv + $repDel + $repIns);
	    my $lengthRD = $rData[4];
	    my $repStartRD = $rData[5];
	    my $repStopRD = $rData[6];
	    my $sizeInterval = $size1{$interval};
	    my $repSizeForID = $repSizes1{$repID};
	    my $repCountsInInterval = $repCountsInterval1{$interval};
	    
	    $info = $repID ."\t". $id ."\t". $id2rep{$id} ."\t". $id2c{$id} ."\t". $id2f{$id}."\t". $interval ."\t". $sizeInterval  ."\t". $repData{$id} ."\t". $pIdentity;
	    $info = $info ."\t". $repSizeForID ."\t". $contigBedRepLen ."\t". $repCountsInInterval ."\t". $totalStats1;

	    if($i == 0){#only print the output file once, still need to loop through to get the other categories of repeats
		$sUsed{$id2rep{$id}} = 1;
		$cUsed{$id2c{$id}} = 1;
		$fUsed{$id2f{$id}} = 1;
		print B1DATA $info ."\n";
	    }

	    #print $info ."\n";
	    my $cat = $id2cat{$id};
	    
	    $b1Intervals{$interval ."--". $whichBed} = 1;
	    my $combo = $interval ."--". $whichBed ."--". $cat;

	    @data = split(/\t/, $info);

	    #pie data collection
	    my $pieKey = $cat ."--". $whichBed;
	    $pieCount{$pieKey}++;
	    $pieLength{$pieKey} += $lengthRD;
	    
	    my $dtLen =  $bed1DNALen + $bed2DNALen;
	    my $rtLen = $bed1RepLen + $bed2RepLen;
	    my $rtCount = $bed1RepCount + $bed2RepCount;
	    my $riCount = $repCountsInInterval;

	    $sumRepLenSeq{$combo} += $lengthRD;
	    $sumRepCountSeq{$combo}++;

	    $sumRepLenNormRepTotalLen{$combo} += $lengthRD / $rtLen;
	    $sumRepLenNormDNATotalLen{$combo} += $lengthRD / $dtLen;
	    $sumRepLenNormRepIntervalLen{$combo} += $lengthRD / $contigBedRepLen;
	    $sumRepLenNormDNAIntervalLen{$combo} += $lengthRD / $sizeInterval;

	    $sumRepCountNormRepTotalLen{$combo} += 1 / $rtLen;
	    $sumRepCountNormDNATotalLen{$combo} += 1 / $dtLen;
	    $sumRepCountNormRepIntervalLen{$combo} += 1 / $contigBedRepLen;
	    $sumRepCountNormDNAIntervalLen{$combo} += 1 / $sizeInterval;

	    $sumRepCountNormRepTotalCount{$combo} += 1 / $rtCount; 
	    #$sumRepCountNormDNATotalCount{$combo} +=;
	    $sumRepCountNormRepIntervalCount{$combo} += 1 / $riCount;
	    #$sumRepCountNormDNAIntervalCount{$combo} += 1 / $data[];

	    #this calculates a "rolling" average
	    if($avgDivRepSeq{$combo}){
		$avgDivRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgDivRepSeq{$combo}) + $repDiv) / $sumRepCountSeq{$combo};
		$avgDelRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgDelRepSeq{$combo}) + $repDel) / $sumRepCountSeq{$combo};
		$avgInsRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgInsRepSeq{$combo}) + $repIns) / $sumRepCountSeq{$combo};
		$avgLenRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgLenRepSeq{$combo}) + $lengthRD) / $sumRepCountSeq{$combo};
		$avgIdtRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgIdtRepSeq{$combo}) + $pIdentity) / $sumRepCountSeq{$combo};
		$avgStopRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgStopRepSeq{$combo}) + $repStopRD) / $sumRepCountSeq{$combo};
		$avgStartRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgStartRepSeq{$combo}) + $repStartRD) / $sumRepCountSeq{$combo};
		#(size * average + value) / (size + 1);
	    }else{
		$avgDivRepSeq{$combo} = $repDiv / $sumRepCountSeq{$combo};
		$avgDelRepSeq{$combo} = $repDel / $sumRepCountSeq{$combo};
		$avgInsRepSeq{$combo} = $repIns / $sumRepCountSeq{$combo};
		$avgLenRepSeq{$combo} = $lengthRD / $sumRepCountSeq{$combo};
		$avgIdtRepSeq{$combo} = $pIdentity / $sumRepCountSeq{$combo};
		$avgStopRepSeq{$combo} = $repStopRD / $sumRepCountSeq{$combo};
                $avgStartRepSeq{$combo} = $repStartRD / $sumRepCountSeq{$combo};

	    }
	}

	#print "on bed2 regions\n";
	my $totalStats2 = $bed2DNALen ."\t". $bed1DNALen ."\t". $bed2RepLen ."\t". $bed1RepLen ."\t". $bed2RepCount ."\t". $bed1RepCount;
	foreach my $repID (@bed2ID){
	    my @data = split(/_on_/, $repID);
	    my $interval = $data[1];
	    my $id = $data[0];
	    my $contigBedRepLen = $b2IntervalRepLen{$interval};

	    my $info;
            my $whichBed = "bed2";
            my @rData = split(/\t/, $repData{$id});
            my $repDiv = $rData[1];
            my $repDel = $rData[2];
            my $repIns = $rData[3];
            my $pIdentity = 100 - ($repDiv + $repDel + $repIns);
            my $lengthRD = $rData[4];
            my $repStartRD = $rData[5];
            my $repStopRD = $rData[6];
            my $sizeInterval = $size2{$interval};
            my $repSizeForID = $repSizes2{$repID};
            my $repCountsInInterval = $repCountsInterval2{$interval};

            $info = $repID ."\t". $id ."\t". $id2rep{$id} ."\t". $id2c{$id} ."\t". $id2f{$id}."\t". $interval ."\t". $sizeInterval  ."\t". $repData{$id} ."\t". $pIdentity;
            $info = $info ."\t". $repSizeForID ."\t". $contigBedRepLen ."\t". $repCountsInInterval ."\t". $totalStats2;

            if($i == 0){#only print the output file once, still need to loop through to get the other categories of repeats
                $sUsed{$id2rep{$id}} = 1;
                $cUsed{$id2c{$id}} = 1;
                $fUsed{$id2f{$id}} = 1;
                print B2DATA $info ."\n";
            }
            #print $info ."\n";
            my $cat = $id2cat{$id};

            $b2Intervals{$interval ."--". $whichBed} = 1;
            my $combo = $interval ."--". $whichBed ."--". $cat;

            @data = split(/\t/, $info);

            #pie data collection
            my $pieKey = $cat ."--". $whichBed;
            $pieCount{$pieKey}++;
            $pieLength{$pieKey} += $lengthRD;

            my $dtLen =  $bed1DNALen + $bed2DNALen;
            my $rtLen = $bed1RepLen + $bed2RepLen;
            my $rtCount = $bed1RepCount + $bed2RepCount;
            my $riCount = $repCountsInInterval;

            $sumRepLenSeq{$combo} += $lengthRD;
            $sumRepCountSeq{$combo}++;

            $sumRepLenNormRepTotalLen{$combo} += $lengthRD / $rtLen;
            $sumRepLenNormDNATotalLen{$combo} += $lengthRD / $dtLen;
            $sumRepLenNormRepIntervalLen{$combo} += $lengthRD / $contigBedRepLen;
            $sumRepLenNormDNAIntervalLen{$combo} += $lengthRD / $sizeInterval;

            $sumRepCountNormRepTotalLen{$combo} += 1 / $rtLen;
            $sumRepCountNormDNATotalLen{$combo} += 1 / $dtLen;
            $sumRepCountNormRepIntervalLen{$combo} += 1 / $contigBedRepLen;
            $sumRepCountNormDNAIntervalLen{$combo} += 1 / $sizeInterval;

            $sumRepCountNormRepTotalCount{$combo} += 1 / $rtCount;
            #$sumRepCountNormDNATotalCount{$combo} +=;
            $sumRepCountNormRepIntervalCount{$combo} += 1 / $riCount;
            #$sumRepCountNormDNAIntervalCount{$combo} += 1 / $data[];

            #this calculates a "rolling" average
            if($avgDivRepSeq{$combo}){
                $avgDivRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgDivRepSeq{$combo}) + $repDiv) / $sumRepCountSeq{$combo};
                $avgDelRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgDelRepSeq{$combo}) + $repDel) / $sumRepCountSeq{$combo};
                $avgInsRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgInsRepSeq{$combo}) + $repIns) / $sumRepCountSeq{$combo};
                $avgLenRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgLenRepSeq{$combo}) + $lengthRD) / $sumRepCountSeq{$combo};
                $avgIdtRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgIdtRepSeq{$combo}) + $pIdentity) / $sumRepCountSeq{$combo};
                $avgStopRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgStopRepSeq{$combo}) + $repStopRD) / $sumRepCountSeq{$combo};
                $avgStartRepSeq{$combo} = ((($sumRepCountSeq{$combo} - 1) * $avgStartRepSeq{$combo}) + $repStartRD) / $sumRepCountSeq{$combo};
                #(size * average + value) / (size + 1);
            }else{
                $avgDivRepSeq{$combo} = $repDiv / $sumRepCountSeq{$combo};
                $avgDelRepSeq{$combo} = $repDel / $sumRepCountSeq{$combo};
                $avgInsRepSeq{$combo} = $repIns / $sumRepCountSeq{$combo};
                $avgLenRepSeq{$combo} = $lengthRD / $sumRepCountSeq{$combo};
                $avgIdtRepSeq{$combo} = $pIdentity / $sumRepCountSeq{$combo};
                $avgStopRepSeq{$combo} = $repStopRD / $sumRepCountSeq{$combo};
                $avgStartRepSeq{$combo} = $repStartRD / $sumRepCountSeq{$combo};

            }
	}
    }
    close(B1DATA);
    close(B2DATA);

    #print lists of repeats, classes and families
    print "writing list of all reapeats to files\n";
    my @sL = sort(keys(%sUsed));
    my @cL = sort(keys(%cUsed));
    my @fL = sort(keys(%fUsed));

    my $sjoin = join("\n", @sL);
    my $cjoin = join("\n", @cL);
    my $fjoin = join("\n", @fL);
    my @repListOut = ($sjoin, $cjoin, $fjoin);
    for(my $i = 0; $i < @repListOut; $i++){
	my $repListTSV = $wd ."/". $categories[$i] ."-repeatList.tsv";
	open(RLIST, ">", $repListTSV);
	print RLIST $repListOut[$i] ."\n";
	close(RLIST);
    }

    #print lists of intervals
    my @b1i = keys(%b1Intervals);
    my @b2i = keys(%b2Intervals);

    my @intervals = (@b1i, @b2i);
    my $iOut = join("\n", @intervals);
    my $iOutFile = $wd . "/intervalLabels-c.tsv";
    open(IOUT, ">", $iOutFile);

    print "writing the matrices for R to read in\n";

    my @statsMatData = (\%avgDivRepSeq, \%avgDelRepSeq, \%avgInsRepSeq, \%avgIdtRepSeq, \%avgLenRepSeq, \%sumRepLenSeq, \%sumRepCountSeq, 
			\%sumRepLenNormRepTotalLen, \%sumRepLenNormDNATotalLen, \%sumRepLenNormRepIntervalLen, \%sumRepLenNormDNAIntervalLen,
			\%sumRepCountNormRepTotalLen, \%sumRepCountNormDNATotalLen, \%sumRepCountNormRepIntervalLen, \%sumRepCountNormDNAIntervalLen,
			\%sumRepCountNormRepTotalCount, \%sumRepCountNormRepIntervalCount, \%avgStartRepSeq, \%avgStopRepSeq );
    
    my @statsList = ("avgDivRepSeq", "avgDelRepSeq", "avgInsRepSeq", "avgIdtRepSeq", "avgLenRepSeq", "sumRepLenSeq", "sumRepCountSeq",
		     "sumRepLenNormRepTotalLen", "sumRepLenNormDNATotalLen", "sumRepLenNormRepIntervalLen", "sumRepLenNormDNAIntervalLen",
		     "sumRepCountNormRepTotalLen", "sumRepCountNormDNATotalLen", "sumRepCountNormRepIntervalLen", "sumRepCountNormDNAIntervalLen",
		     "sumRepCountNormRepTotalCount", "sumRepCountNormRepIntervalCount", "avgStartRepSeq", "avgStopRepSeq");
    

    my $printB1 = join("\t", @b1i);
    my $printB2 = join("\t", @b2i);
    print IOUT $printB1 ."\t". $printB2 ."\n";
    close IOUT;

    #this is where the ROI vs background averages can be calculated and printed to a file
    my @repclassList = (\@sL, \@cL, \@fL);
    my @classLabel = ("species", "class", "family");

    print "writing data for bed1\n";
    for(my $j = 0; $j < @statsList; $j++){
	my $stat = $statsList[$j];
	my %matData = %{$statsMatData[$j]};
	#print "writing data for $stat in bed1\n";
	for(my $i = 0; $i < @classLabel; $i++){
	    my $b1MatFile = $wd ."/". $stat ."-". $classLabel[$i] ."-b1.tsv";
	    open(B1MAT, ">", $b1MatFile);
	    #print $b1MatFile ."\n";
	    foreach my $reg (@b1i){
		my @data;
		foreach my $r (@{$repclassList[$i]}){
		    my $combo = $reg ."--". $r;
		    #print $combo ."\n";
		    if(1){
			my $v = "NaN";
			if($matData{$combo}){
			    $v = $matData{$combo};
			}
			push(@data, $v);
		    }else{
			if($matData{$combo}){
			    my $v = $matData{$combo};
			    push(@data, $v);
			}
		    }
		}
		my $line = join("\t", @data);
		print B1MAT $line ."\n";
	    }
	    close(B1MAT);
	}
    }

    print "writing data for bed1\n";
    for(my $j = 0; $j < @statsList; $j++){
	my $stat = $statsList[$j];
	my %matData = %{$statsMatData[$j]};
	#print "writing data for $stat in bed2\n";
	for(my $i = 0; $i < @classLabel; $i++){
	    my $b2MatFile = $wd ."/". $stat ."-". $classLabel[$i] ."-b2.tsv";
	    open(B2MAT, ">", $b2MatFile);
	    #print $b2MatFile ."\n";
	    foreach my $reg (@b2i){
		my @data;
		foreach my $r (@{$repclassList[$i]}){
		    my $combo = $reg ."--". $r;
		    #print $combo ."\n";
		    if(1){
			my $v = "NaN";
			if($matData{$combo}){
			    $v = $matData{$combo};
			}
			push(@data, $v);
		    }else{
			if($matData{$combo}){
			    my $v = $matData{$combo};
			    push(@data, $v);
			}
		    }
		}
		
		my $line = join("\t", @data);
		print B2MAT $line ."\n";
	    }
	    close(B2MAT);
	}
	
    }
    
    print "printing summary data files for pie charts\n";
    my %pieCountNorm2DNA;
    
    my @pieList = ("countRaw", "bpRaw", "countNorm2Rep", "bpNorm2Rep", "countNorm2DNA", "bpNorm2DNA");
    
    #this is where the ROI vs background averages can be calculated and printed to a file

    #pie data collection
    #my $pieKey = $cat ."--bed1";
    #$pieCount{$pieKey}++;
    #$pieLength{$pieKey} += $data[4];
    
    #repeats lists
    #my @sL = sort(keys(%sUsed));
    #my @cL = sort(keys(%cUsed));
    #my @fL = sort(keys(%fUsed));
    
    
    my @repPieList = (\@sL, \@fL, \@cL);
    
    for(my $i = 0; $i < @classLabel; $i++){
	for(my $j = 0; $j < 6; $j = $j + 2){
	    my $statCount = $pieList[$j];
	    my $statBP = $pieList[$j + 1];
	    
	    open( PIECOUNT, ">", $wd ."/pie-". $classLabel[$i] ."-". $statCount .".tsv");
	    open( PIEBP, ">", $wd ."/pie-". $classLabel[$i] ."-". $statBP .".tsv");
	    
	    print "writing ". $classLabel[$i] ." for ". $pieList[$j] ." and $pieList[$j+1]\n";
	    
	    my @rpl = @{$repPieList[$i]};
	    foreach my $pieRep (@rpl){
		my @bedTypes = ("--bed1", "--bed2");
		foreach my $bt (@bedTypes){
		    if(!$pieCount{$pieRep . $bt}){
			$pieCount{$pieRep . $bt} = 0;
		    }
		    if(!$pieLength{$pieRep . $bt}){
			$pieLength{$pieRep . $bt} = 0;
		    }
		}
		
		
		if($j == 0){
		    print PIECOUNT $pieRep ."\t". $pieCount{$pieRep ."--bed1"} ."\t". $pieCount{$pieRep ."--bed2"} ."\n";
		    print PIEBP $pieRep ."\t". $pieLength{$pieRep ."--bed1"} ."\t". $pieLength{$pieRep ."--bed2"} ."\n";
		}elsif($j == 2){ #normalized to rep length
		    print PIECOUNT $pieRep ."\t". ( $pieCount{$pieRep ."--bed1"} / $bed1RepLen ) ."\t". ( $pieCount{$pieRep ."--bed2"} / $bed2RepLen ) ."\n";
		    print PIEBP $pieRep ."\t". ( $pieLength{$pieRep ."--bed1"} / $bed1RepLen ) ."\t". ( $pieLength{$pieRep ."--bed2"} / $bed2RepLen ) ."\n";
		}elsif($j == 4){ #normalized to dna length
		    print PIECOUNT $pieRep ."\t". ( $pieCount{$pieRep ."--bed1"} / $bed1DNALen ) ."\t". ( $pieCount{$pieRep ."--bed2"} / $bed2DNALen ) ."\n";
		    print PIEBP $pieRep ."\t". ( $pieLength{$pieRep ."--bed1"} / $bed1DNALen ) ."\t". ( $pieLength{$pieRep ."--bed2"} / $bed2DNALen ) ."\n";
		}
	    }
	    
	    if($j == 0){
		print PIEBP "Non-repetitive\t". ($bed1DNALen - $bed1RepLen) ."\t". ($bed2DNALen - $bed2RepLen) ."\n";
	    }elsif($j == 2){ #normalized to rep length
		#skip it because it doesn't make sense
	    }elsif($j == 4){ #normalized to dna length
		print PIEBP "Non-repetitive\t". (($bed1DNALen - $bed1RepLen) / $bed1DNALen) ."\t". (($bed2DNALen - $bed2RepLen) / $bed2DNALen) ."\n";
	    }
	    
	    close(PIECOUNT);
	    close(PIEBP);
	}
    }
    
    my $matStatFile = $wd . "/matStatsList.tsv";
    open(MSOUT, ">", $matStatFile);
    foreach my $matStat (@statsList){
	print MSOUT $matStat ."\n";
    }
}
my $rFile = $wd . "/stats.R";
open(ROUT, ">", $rFile);

print ROUT '

dirName <- c("'. $wd .'")
nCores <- '. $cores .'
goIndividual <- '. $runITests .'
goCluster <- '. $runClustTests .'
goClassify <- '. $runClassTests .'
bed1Name <- c("'. $filename1 .'")
bed2Name <- c("'. $filename2 .'")
lengthMultiply <- c("'. $multiLength .'");
iMetrics <- c("'. $runIMetrics .'")
runGT <- '. $runGTTests .'
runLT <- '. $runLTTests .'
run2S <- '. $run2STests .'


tailTypes <- c(run2S, runGT, runLT)

if(0){
    nCores <- 20
    goIndividual <- 1
    goCluster <- 1
    goClassify <- 1
    bed1Name <- c("SRR5723785.fastq.fa.15k")
    bed2Name <- c("SRR5723786.fastq.fa.15k")
    mcLength <- c("100000");
    iMetrics <- c("T,T,T,T,T,T,T,T,T,T,T")
    tailTypes <- c(0,0,1)
}

setwd(dirName)

#writeDir <- paste(tempdir(), dirName, sep="/")
writeDir <-  paste("working", dirName, sep="/")
writeDir <-  "."
dir.create(writeDir, showWarnings = FALSE)

iMetrics <- as.logical(unlist(strsplit(iMetrics, split=",")))
iTests <- 1:11

options(max.print=1000000)

library(ggbiplot)
library(ggplot2)
library(Matching)
library(kSamples)
library(parallel)
library(gplots)
library(ape)
library(data.table)

library(magrittr)
library(vegan)
library(pvclust)
library(mclust)
library(fpc)
library(class)
library(Hmisc)
library(coin)

###following variables usesful for debugging
b1File <- c("bed1Stats.tsv")
b2File <- c("bed2Stats.tsv")
###read in data
b1FullData <- as.data.frame(fread(b1File, header=FALSE, sep = "\t"))
b2FullData <- as.data.frame(fread(b2File, header=FALSE, sep = "\t"))

catLabels <- c("species", "class", "family")
repeatFileList <- c("species-repeatList.tsv", "class-repeatList.tsv", "family-repeatList.tsv")

regionList <- as.data.frame(fread("intervalLabels-c.tsv", header=FALSE, sep = "\t"))
regionList <- sapply(regionList, as.character)

###get the regions that are overlapped, (search for the ::: in contig:::start-stop)
rmRG <- grepl ("--bed1", regionList)
bed1Reg <- regionList[rmRG]
bed2Reg <- regionList[!rmRG]

#options(warn=1)

createFigures <- function(statsTest, repLevel, repName, rmb1, rmb2, writeDir){
    b1CDF <- ecdf(rmb1)
    b2CDF <- ecdf(rmb2)

    compressSwitch <- FALSE

    figFileName <- paste(writeDir, "/individualRepeatTests/", repLevel,"/", repLevel, "-", statsTest,"/cdf/CDF-", statsTest, "-for-", repName, ".png", sep = "")
    png(file = figFileName, width = 1200, height = 1200, units = "px", type = "cairo")
    #pdf(file = paste(writeDir, "/individualRepeatTests/", repLevel,"/", repLevel, "-", statsTest,"/cdf/CDF-", statsTest, "-for-", repName, ".pdf", sep = ""), compress=compressSwitch)

    plot(b1CDF, col="blue", main= paste(statsTest, " for ", repName, " (", repLevel, ") ", sep = ""), xlab=statsTest)
    lines(b2CDF, col="red")
    legend(x="bottomright", legend= paste("blue =", bed1Name, " regions\nred =", bed2Name , " regions", sep = " "))
    dev.off()

    b1df <- as.data.frame(rmb1)
    b1df$type <- bed1Name
    colnames(b1df) <- c("stat", "bed")
    b2df <- as.data.frame(rmb2)
    b2df$type <- bed2Name
    colnames(b2df) <- c("stat", "bed")
    bcdf <- rbind(b1df, b2df)
    
    figFileName <- paste(writeDir, "/individualRepeatTests/", repLevel,"/", repLevel, "-", statsTest,"/violin/violin-", statsTest, "-for-", repName, ".png", sep = "")
    if((length(rmb1) > 2) && (length(rmb2) > 2)){
        #pdf(file = figFileName, compress=compressSwitch)
        png(file = figFileName, width = 1200, height = 1200, units = "px", type = "cairo")
        p <- ggplot(bcdf, aes(x=bed, y=stat, color=bed)) + ylab(statsTest) + geom_violin()
        p <- p + stat_summary(fun.data="mean_sdl", geom="pointrange")
        p <- p + scale_color_manual(values=c("#0000FF", "#FF0000"))
        print(p)
        dev.off()
    }else if ((length(rmb1) > 2) && (length(rmb2) < 3)){
        #pdf(file = figFileName, compress=compressSwitch)
        png(file = figFileName, width = 1200, height = 1200, units = "px", type = "cairo")
        p <- ggplot(b1df, aes(x=bed, y=stat, color=bed)) + ylab(statsTest) + geom_violin()
        p <- p + stat_summary(fun.data="mean_sdl", geom="pointrange")
        p <- p + scale_color_manual(values=c("#FF0000"))
        print(p)
        dev.off()
    }else if ((length(rmb1) < 3) && (length(rmb2) > 2)){
        #pdf(file = figFileName, compress=compressSwitch)
        png(file = figFileName, width = 1200, height = 1200, units = "px", type = "cairo")
        p <- ggplot(b2df, aes(x=bed, y=stat, color=bed)) + ylab(statsTest) + geom_violin()
        p <- p + stat_summary(fun.data="mean_sdl", geom="pointrange")
        p <- p + scale_color_manual(values=c("#0000FF"))
        print(p)
        dev.off()
    }
}

onlyQuickTests <- function(curCombo, allCombos, statsTypes, statsNames, catLabels, tailTypes, b1, b2, writeDir){
    #curCombo <- 14873
    combo <- allCombos[curCombo,]
    
    stat <- as.integer(paste(unlist(combo[2]), collapse=""))
    repName <- paste(unlist(combo[1]), collapse="")
    repLevelNum <- as.integer(combo[3])
    repLevel <- catLabels[repLevelNum]

    rmb1 <- b1[which(b1$repeatLevel == repLevelNum & b1$repeatName == repName),stat+1]
    rmb2 <- b2[which(b2$repeatLevel == repLevelNum & b2$repeatName == repName),stat+1]

    b1Len <- length(rmb1)
    b2Len <- length(rmb2)

    repTestPrintInfo <- paste(repName, "--", statsTypes[stat], "with bed1 count:", length(rmb1),"- bed2 count:",  length(rmb2), sep=" ")

    if(curCombo %% 1000 == 0){
        ### write("prints to stderr", stderr())
	outputLine <- paste(Sys.time(), "-- phase 1/2 repeat test:", curCombo, "/", dim(allCombos)[1], "-- test of:", repTestPrintInfo, sep=" ")
        print(outputLine)
        ### stderr(outputLine)
    }

    nPTests <- 7 #- sum(tailTypes) - sum(tailTypes)

    rmDF <- data.frame(vals = c(rmb1, rmb2), bed = rep(c("bed1","bed2"), c(length(rmb1), length(rmb2))))

    #one repeat is not found in sample
    if((length(rmb1) == 0) || (length(rmb2) == 0)){
        return(c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,rep(0, nPTests)))
    }

    if(length(unique(rmb1)) == 1 && length(unique(rmb2)) == 1){
        if(unique(rmb1) == unique(rmb2)){
            return(c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,rep(0, nPTests)))
        }
    }

    ### other possibilities Cramer-von Mises, and ryan-joiner

    if((b1Len + b2Len) > 500 ) {
        #print(paste(Sys.time(), "-- skipping MWU/WSR/AD tests:", curCombo, "/", dim(allCombos)[1], "-- test type:", statsTypes[stat], "--", repName,  "with bed1 count:", length(rmb1),"- bed2 count:",  length(rmb2), sep=" "))
        return(c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,rep(-1, nPTests)))
    }

    adRes <- ad.test(rmb1, rmb2)
    ad.pvalue <- adRes$ad[2,3]
   
    if((ad.pvalue < 1 * (10 ^ -5)) | ((1 - ad.pvalue) < 1 * (10 ^ -5))){ #ad.test is only accurate within a certain range of pvalues
        #print(paste(Sys.time(), "-- skipping MWU/WSR/AD tests:", curCombo, "/", dim(allCombos)[1], "-- test type:", statsTypes[stat], "--", repName,  "with bed1 count:", length(rmb1),"- bed2 count:",  length(rmb2), sep=" "))
        return(c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,rep(-1, nPTests)))
    }

    numEstimateSimulationsWT <- 1 * (10 ^ 5) + 1 
    numSimulationsWorstCaseWT <- 5 * (10 ^ 7) + 1
    numSimulationsBestCaseWT <- 5 * (10 ^ 5) + 1

    mwu.t <- -2
    mwu.g <- -2
    mwu.l <- -2

    mwu.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsBestCaseWT)))[1]
    
    if((mwu.t < 1/(numSimulationsBestCaseWT * 20)) | ((1 - mwu.t) < 1/(numSimulationsBestCaseWT * 20))){
        ### the pvalue is so extreme that running in parallel for each repeat will reach an accurate p-value
        return(c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,rep(-1, nPTests)))
    }else{
        numSimulationsMWU <- numSimulationsBestCaseWT
        if(tailTypes[2]){
            mwu.g <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU), alternative="greater"))[1]
        }
        if(tailTypes[3]){ 
            mwu.l <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU), alternative="less"))[1]
        }
    }
    if(!(tailTypes[1])){
        mwu.t <- -2
    }
    mwu.pvalues <- c(mwu.t, mwu.g, mwu.l)
    
    wsr.t <- -2
    wsr.g <- -2
    wsr.l <- -2

    wsr.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsBestCaseWT)))[1]
    if((wsr.t < 1/(numSimulationsBestCaseWT * 20)) | ((1 - wsr.t) < 1/(numSimulationsBestCaseWT * 20))){
        return(c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,rep(-1, nPTests)))
    }else{
        numSimulationsWSR <- numSimulationsBestCaseWT
        if(tailTypes[2]){ 
            wsr.g <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR), paired=TRUE, alternative="greater"))[1]
        }
        if(tailTypes[3]){ 
            wsr.l <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR), paired=TRUE, alternative="less"))[1]
        }
    }
    if(!(tailTypes[1])){
        wsr.t <- -2
    }

    wsr.pvalues <- c(wsr.t, wsr.g, wsr.l)

    res <- c(ad.pvalue, mwu.pvalues, wsr.pvalues)
    
    smallSelection <- which(res >= 0 & res < 0.1)
    largeSelection <- which(res > 0.9 & res <= 1)
    
    if(length(c(smallSelection, largeSelection)) > 0){
        createFigures(statsTest = statsTypes[stat], repLevel = repLevel, repName = repName, rmb1 = rmb1, rmb2 = rmb2, writeDir = writeDir)
        #print(paste(Sys.time(), "-- creating figure for:", curCombo, "/", dim(allCombos)[1], "-- test type:", statsTypes[stat], "--", repName,  "with bed1 count:", length(rmb1),"- bed2 count:",  length(rmb2), sep=" "))
    }
    #print(paste(Sys.time(), "-- finished :", curCombo, "/", dim(allCombos)[1], "-- test of:", repTestPrintInfo, sep=" "))
    res[res==0] <- 2.2e-16
    results <- c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,res)
    return(results)
}

adTestMulticore <- function(thread, numIterations, rmb1, rmb2){
    adRes <- ad.test(rmb1, rmb2, method= "simulated", Nsim = numIterations, dist=FALSE)
    ad.pvalue <- adRes$ad[2,4]
    numObserved <- numIterations * ad.pvalue
    return(numObserved)
}

onlyLongTests <- function(longCombos, allCombos, statsTypes, statsNames, catLabels, tailTypes, b1, b2, writeDir, pThreads){
    #curCombo <- 14873
    nOtherCol <- 6
    nPTests <- 7
    longResults <- matrix(0, nrow = length(longCombos), ncol = (nOtherCol + nPTests))

    for(curLong in 1:length(longCombos)){
        curCombo <- as.integer(longCombos[curLong])
        combo <- allCombos[curCombo,]
        #combo <- partialAllCombos[curCombo,]

        stat <- as.integer(paste(unlist(combo[2]), collapse=""))
        repName <- paste(unlist(combo[1]), collapse="")
        repLevelNum <- as.integer(combo[3])
        repLevel <- catLabels[repLevelNum]

        rmb1 <- b1[which(b1$repeatLevel == repLevelNum & b1$repeatName == repName),stat+1]
        rmb2 <- b2[which(b2$repeatLevel == repLevelNum & b2$repeatName == repName),stat+1]
        b1Len <- length(rmb1)
        b2Len <- length(rmb2)
        rmDF <- data.frame(vals = c(rmb1, rmb2), bed = rep(c("bed1","bed2"), c(length(rmb1), length(rmb2))))
        expBC <- 6

        repTestPrintInfo <- paste("Testing", repName, "for", statsNames[stat], "", stat, "/", length(statsNames), ") with counts:", length(rmb1),"and",  length(rmb2), sep=" ")

        if(curLong %% 5 == 0){
            outputLine <- paste(Sys.time(), "-- phase 2/2 tests:", curLong, "/", length(longCombos), "--", repTestPrintInfo, sep=" ")
            print(outputLine)
            ### stderr(outputLine)
        }

        ###consider looking into Cramer-von Mises, Wilcoxon-Mann-Whitney and ryan-joiner
        adRes <- ad.test(rmb1, rmb2)
        ad.pvalue <- adRes$ad[2,3]

        if((ad.pvalue < 1 * (10 ^ -5)) | ((1 - ad.pvalue) < 1 * (10 ^ -5))){ #ad.test is only accurate within a certain range of pvalues
            numIterations <- 1 * (10 ^ expBC)

            adEstimateSims <- 5000
            numSubEstimateIterations <- ceiling(adEstimateSims / pThreads)

            simStartTime <- as.numeric(as.POSIXct(Sys.time()))
            numObserved <- mclapply(1:pThreads, adTestMulticore, numIterations = numSubEstimateIterations, rmb1=rmb1, rmb2=rmb2, mc.cores=pThreads, mc.preschedule=TRUE)
            simStopTime <- as.numeric(as.POSIXct(Sys.time()))
            simSpentTime <- simStopTime - simStartTime
            estimatedTimeMinutesAD <- signif(((simSpentTime * (numIterations / adEstimateSims)) / 60), 3)
            if(estimatedTimeMinutesAD > 10){
                outputLine <- paste(Sys.time(), "-- Long Calculation Notification - AD simulation test:", curLong, "/", length(longCombos), "-- estimating", estimatedTimeMinutesAD, "minutes -- test of:", repTestPrintInfo, sep=" ")
                print(outputLine)
                ### stderr(outputLine)

            }
            numSubIterations <- ceiling(numIterations / pThreads)
            simStartRealTime <- as.numeric(as.POSIXct(Sys.time()))
            numObserved <- mclapply(1:pThreads, adTestMulticore, numIterations = numSubIterations, rmb1=rmb1, rmb2=rmb2, mc.cores=pThreads, mc.preschedule=TRUE)
            simStopRealTime <- as.numeric(as.POSIXct(Sys.time()))
            simSpentRealTime <- signif(((simStopRealTime - simStartRealTime)/60), 3)

            countObserved <- sum(unlist(numObserved))
            ad.pvalue <- countObserved / numIterations
            if(ad.pvalue == 0){
                ad.pvalue <- 1/(numIterations + 1)
            }
            if(ad.pvalue == 1){
                ad.pvalue <- numIterations/(numIterations + 1)
            }
        }

        numEstimateSimulationsWT <- 1 * (10 ^ 3) + 1 
        numSimulationsWorstCaseWT <- 1 * (10 ^ 7) + 1
        numSimulationsBestCaseWT <- 1 * (10 ^ 5) + 1

        mwu.t <- -2
        mwu.g <- -2
        mwu.l <- -2

        #print(paste(Sys.time(), "-- phase 2/2 MWU tests:", curLong, "/", length(longCombos), sep=" "))

        simStartTimeMWU <- as.numeric(as.POSIXct(Sys.time()))
        mwu.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsBestCaseWT, parallel = "multicore", ncpus = pThreads)))[1]
        simStopTimeMWU <- as.numeric(as.POSIXct(Sys.time()))
        simSpentTimeMWU <- simStopTimeMWU - simStartTimeMWU

        if((mwu.t < 1/(numSimulationsBestCaseWT * 20)) | ((1 - mwu.t) < 1/(numSimulationsBestCaseWT * 20))){
            numSimulationsMWU <- numSimulationsWorstCaseWT
            estimatedTimeMinutesMWU <- signif(((sum(tailTypes) * (simSpentTimeMWU * (numSimulationsMWU / numSimulationsBestCaseWT)) / 60)), 3)
            if(estimatedTimeMinutesMWU > 10){
                outputLine <- paste(Sys.time(), "-- Long Calculation Notification - MWU:", curLong, "/", length(longCombos), "-- estimating", estimatedTimeMinutesMWU, "minutes -- test of:", repTestPrintInfo, sep=" ")
                print(outputLine)
                ### stderr(outputLine)
            }
            if(tailTypes[1]){
                mwu.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU, parallel = "multicore", ncpus = pThreads)))[1]
            }
            if(tailTypes[2]){
                mwu.g <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU, parallel = "multicore", ncpus = pThreads), alternative="greater"))[1]
            }
            if(tailTypes[3]){
                mwu.l <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU, parallel = "multicore", ncpus = pThreads), alternative="less"))[1]
            }
        }else{
            numSimulationsMWU <- numSimulationsBestCaseWT
            estimatedTimeMinutesMWU <- signif((((tailTypes[2] + tailTypes[3]) * (simSpentTimeMWU * (numSimulationsMWU / numSimulationsBestCaseWT)) / 60)), 3)
            if(estimatedTimeMinutesMWU > 10){
                outputLine <- paste(Sys.time(), "-- Long Calculation Notification - MWU:", curLong, "/", length(longCombos), "-- estimating", estimatedTimeMinutesMWU, "minutes -- test of:", repTestPrintInfo, sep=" ")
                print(outputLine)
                ### stderr(outputLine)
            }
            #mwu.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU, parallel = "multicore", ncpus = pThreads)))[1]
            if(tailTypes[2]){
                mwu.g <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU, parallel = "multicore", ncpus = pThreads), alternative="greater"))[1]
            }
            if(tailTypes[3]){
                mwu.l <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsMWU, parallel = "multicore", ncpus = pThreads), alternative="less"))[1]
            }
        }
        if(!(tailTypes[1])){
            mwu.t <- -2
        } 
        mwu.pvalues <- c(mwu.t, mwu.g, mwu.l)

        wsr.t <- -2
        wsr.g <- -2
        wsr.l <- -2

        #print(paste(Sys.time(), "-- phase 2/2 WSR tests:", curLong, "/", length(longCombos), sep=" "))
        simStartTimeWSR <- as.numeric(as.POSIXct(Sys.time()))
        wsr.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsBestCaseWT, parallel = "multicore", ncpus = pThreads)))[1]
        simStopTimeWSR <- as.numeric(as.POSIXct(Sys.time()))
        simSpentTimeWSR <- simStopTimeWSR - simStartTimeWSR
        
        if((wsr.t < 1/(numSimulationsBestCaseWT * 20)) | ((1 - wsr.t) < 1/(numSimulationsBestCaseWT * 20))){
            numSimulationsWSR <- numSimulationsWorstCaseWT
            estimatedTimeMinutesWSR <- signif(((sum(tailTypes) * (simSpentTimeWSR * (numSimulationsWSR / numSimulationsBestCaseWT)) / 60)), 3)
            if(estimatedTimeMinutesWSR > 10){
                outputLine <- paste(Sys.time(), "-- Long Calculation Notification - WSR:", curLong, "/", length(longCombos), "-- estimating", estimatedTimeMinutesWSR, "minutes -- test of:", repTestPrintInfo, sep=" ")
                print(outputLine)
                ### stderr(outputLine)
            }
            if(tailTypes[1]){
                wsr.t <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR, parallel = "multicore", ncpus = pThreads), paired=TRUE))[1]
            }
            if(tailTypes[2]){
                wsr.g <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR, parallel = "multicore", ncpus = pThreads), paired=TRUE, alternative="greater"))[1]
            }
            if(tailTypes[3]){
                wsr.l <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR, parallel = "multicore", ncpus = pThreads), paired=TRUE, alternative="less"))[1]
            }

        }else{
            numSimulationsWSR <- numSimulationsBestCaseWT
            estimatedTimeMinutesWSR <- signif((((tailTypes[2] + tailTypes[3]) * (simSpentTimeWSR * (numSimulationsWSR / numSimulationsBestCaseWT)) / 60)), 3)
            if(estimatedTimeMinutesWSR > 10){
                outputLine <- paste(Sys.time(), "-- Long Calculation Notification - WSR:", curLong, "/", length(longCombos), "-- estimating", estimatedTimeMinutesWSR, "minutes -- test of:", repTestPrintInfo, sep=" ")
                print(outputLine)
                ### stderr(outputLine)
            }
            if(tailTypes[2]){
                wsr.g <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR, parallel = "multicore", ncpus = pThreads), paired=TRUE, alternative="greater"))[1]
            }
            if(tailTypes[3]){
                wsr.l <- pvalue(wilcox_test(vals~bed, data=rmDF, distribution = approximate(B = numSimulationsWSR, parallel = "multicore", ncpus = pThreads), paired=TRUE, alternative="less"))[1]
            }
        }
        if(!(tailTypes[1])){
            wsr.t <- -2
        }
        wsr.pvalues <- c(wsr.t, wsr.g, wsr.l)

        res <- c(ad.pvalue, mwu.pvalues, wsr.pvalues)
        #print(paste(Sys.time(), "-- phase 2/2 done with tests:", curLong, "/", length(longCombos), sep=" "))

        smallSelection <- which(res >= 0 & res < 0.1)
        largeSelection <- which(res > 0.9 & res <= 1)

        if(length(c(smallSelection, largeSelection)) > 0){
            createFigures(statsTest = statsTypes[stat], repLevel = repLevel, repName = repName, rmb1 = rmb1, rmb2 = rmb2, writeDir = writeDir)
        }
        
        ##colnames(resultStats) <- c("repeat", "ROI counts", "background counts", "pvalues...")
        #res <- c(ks.pvalues, mwu.pvalues, wsr.pvalues, ad.pvalue)
        res[res==0] <- 2.2e-16
        longResults[curLong,] <- c(stat,repLevelNum,repName,b1Len,b2Len,curCombo,res)
    }

    return(longResults)
}

adjustPValues <- function(repeatDF, pColNames){
    #repeatDF <- partialRawResultsDF
    countCutoffs <- c(7, 8, 0, 10, 20) #7 & 8 correspond to the columns with the tenth of max and hundredth of max count each
    adjCutoffNames <- c("tenth", "hundredth", "zero", "minimumOfTen", "minimumOfTwenty")
    repDF <- repeatDF
    for(co in 1:length(countCutoffs)){
        bed1Counts <- as.integer(as.character(repDF[,5]))
        bed2Counts <- as.integer(as.character(repDF[,6]))
        if(co < 3){
            shortStatsData <- repDF[which((bed1Counts >= as.integer(as.character(repDF[,countCutoffs[co]]))) | (bed2Counts >= as.integer(as.character(repDF[,countCutoffs[co]])))),]
        }else{
            cco <- countCutoffs[co]
            shortStatsData <- repDF[which((bed1Counts >= cco) | ( bed2Counts >= cco)), ]
        }

        if(dim(repDF)[1] > 0){
            for(pTestType in 1:length(pColNames)){
                pvToAdj <- as.numeric(as.character(shortStatsData[,(pTestType + 8)]))
                p.holm <- p.adjust(pvToAdj, method="holm")
                p.bh <- p.adjust(pvToAdj, method="BH")
            
                toMerge <- cbind(p.holm, p.bh)

                adjNames <- c("holm.adjP", "BH.adjP")
                comboAdjNames <- unlist(lapply(adjNames, paste, adjCutoffNames[co], sep="-"))
                comboAdjNames <- unlist(lapply(comboAdjNames, paste, pColNames[pTestType], sep="-"))

                colnames(toMerge) <- comboAdjNames
                toMerge <- as.data.frame(toMerge)
                toMerge$unique_name <- row.names(shortStatsData)
                row.names(toMerge) <- row.names(shortStatsData)
                repDF <- merge(repDF, toMerge, by="unique_name", all=TRUE)
                row.names(repDF) <- as.character(repDF$unique_name)
            }
        }
    }
    return(repDF)
}

###do the individual tests for each level of repeat
if(goIndividual){
    statsTypes <- c("percent_divergence", "percent_deleted", "percent_inserted", "length", "bpRepStart", "bpRepStop", "lenRegRepNorm", "lenRegDNANorm", "lenAllRepNorm", "lenAllDNANorm", "percent_identity")
    statsNames <- c("percent divergence", "percent deleted", "percent inserted", "length", "first base of repeat", "last base of repeat", "fraction of repetitive nucleotides in region", "fraction of DNA in region", "fraction of total repetitive nucleotides", "fraction of total nucleotides", "percent identity")

    allCombos <- matrix(, ncol=3) 
    colnames(allCombos) <- c("repeatName", "testType", "repeatLevel")

    bedStatsColNames <- c("repeatName", statsTypes, "repeatLevel")
    b1 <- matrix(, ncol=13) 
    b2 <- matrix(, ncol=13) 
    colnames(b1) <- bedStatsColNames
    colnames(b2) <- bedStatsColNames

    for ( i in 3:1){
        repType <- catLabels[i]

        ###create plots for individual types of repeats
        ###remove empty elements
        reps <- t(as.data.frame(fread(repeatFileList[i], header=FALSE, sep = "\t")))
        repLen <- length(reps)
        repsName <- reps[1:repLen]
        repeatNames <- reps
        ###create CDF plot for each repeat
        dir.create(paste(writeDir, "/individualRepeatTests/", sep=""), showWarnings = FALSE)
        dir.create(paste(writeDir, "/individualRepeatTests/", repType, sep=""), showWarnings = FALSE)

        ###create CDF and scatter plots for divergence 3, del 4, ins 5, len 6, repStart 7, repStop 8
        ###this is the portion of the script that takes the longest and does not need to be rerun after you
        ###identify the repeats that are of interest
        offset <- i + 2

        b1Build <- b1FullData[offset]
        b2Build <- b2FullData[offset]

        for( st in iTests[iMetrics]){###this is for the columns in the data files
            dir.create(paste(writeDir, "/individualRepeatTests/", repType, "/", repType, "-", statsTypes[st], sep = ""), showWarnings = FALSE)
            dir.create(paste(writeDir, "/individualRepeatTests/", repType, "/", repType, "-",  statsTypes[st],"/violin", sep = ""), showWarnings = FALSE)
            dir.create(paste(writeDir, "/individualRepeatTests/", repType, "/", repType, "-",  statsTypes[st],"/cdf", sep = ""), showWarnings = FALSE)
        }

        ##see bed1Stats.tsv for organization of b1FullData
        for(st in 1:3){
            k <- st + 8
            ##get data for each stats type, pay attention to the offset
            b1Build <- cbind(b1Build,b1FullData[k])
            b2Build <- cbind(b2Build,b2FullData[k])
            cdfStatName <- paste(catLabels[i], statsTypes[st], sep = "-")
            humanName <- statsNames[st]
        }
        ##for length
        b1Build <- cbind(b1Build,b1FullData[12])
        b2Build <- cbind(b2Build,b2FullData[12])
        ##for repStart
        b1Build <- cbind(b1Build,b1FullData[13])
        b2Build <- cbind(b2Build,b2FullData[13])
        ##for repStop
        b1Build <- cbind(b1Build,b1FullData[14])
        b2Build <- cbind(b2Build,b2FullData[14])
        ##for normalized to rep length
        b1Build <- cbind(b1Build,(b1FullData[12]/b1FullData[17]))
        b2Build <- cbind(b2Build,(b2FullData[12]/b2FullData[17]))
        ##for normalized to dna length
        b1Build <- cbind(b1Build,(b1FullData[12]/b1FullData[7]))
        b2Build <- cbind(b2Build,(b2FullData[12]/b2FullData[7]))
        ##for normalized to total Rep length
        b1Build <- cbind(b1Build,(b1FullData[12]/b1FullData[21]))
        b2Build <- cbind(b2Build,(b2FullData[12]/b2FullData[21]))
        ##for normalized to total DNA length
        b1Build <- cbind(b1Build,(b1FullData[12]/b1FullData[19]))
        b2Build <- cbind(b2Build,(b2FullData[12]/b2FullData[19]))
        ##for total percent identity
        b1Build <- cbind(b1Build,b1FullData[15])
        b2Build <- cbind(b2Build,b2FullData[15])

        b1Partial <- cbind(b1Build, i)
        b2Partial <- cbind(b2Build, i)

        colnames(b1Partial) <- bedStatsColNames
        colnames(b2Partial) <- bedStatsColNames

        allCombosPartial <- expand.grid(repeatNames, iTests[iMetrics], i)
        colnames(allCombosPartial) <- c("repeatName", "testType", "repeatLevel")

        allCombos <- rbind(allCombos, allCombosPartial)

        b1 <- rbind(b1, b1Partial)
        b2 <- rbind(b2, b2Partial)
    }

    allCombos <- allCombos[-1,]
    b1 <- b1[-1,]
    b2 <- b2[-1,]
    
    #options(warn=1)

    b1colName <- paste(bed1Name, "counts", sep = "_")
    b2colName <- paste(bed2Name, "counts", sep = "_")

    pKSColNames <- c("KS 2-sided p", "KS greater p", "KS lesser p")
    pMWUColNames <- c("MWU 2-sided p", "MWU greater p", "MWU lesser p")
    pWSRColNames <- c("WSR 2-sided p", "WSR greater p", "WSR lesser p")
    pADColNames <- c("AD discontinuous p")
    #pColNames <- c(pKSColNames, pMWUColNames, pWSRColNames, pADColNames)
    pColNames <- c(pADColNames, pMWUColNames, pWSRColNames)
    #pColNames <- c(pADColNames)
    tailSelection <- c(1, tailTypes, tailTypes)
    pColNames <- pColNames[as.logical(tailSelection)]
    
    if(0){
        #allCombos <- expand.grid(repeatNames[1:10], 1:11)
        pThreads <- 30
        mcThreads <- 30
        printThreads <- as.integer(nCores / 2)

        partialAllCombos <- allCombos[seq(from=1, to=dim(allCombos)[1], by=23),printThreads]

        system.time(statsPrints <- mclapply(1:dim(partialAllCombos)[1], runCreateFigures, allCombos = partialAllCombos, statsTypes = statsTypes, statsNames = statsNames, catLabels = catLabels, b1=b1, b2=b2, mc.cores=printThreads, mc.preschedule=FALSE))
        system.time(statsDataResQuick <- mclapply(1:dim(partialAllCombos)[1], onlyQuickTests, allCombos = partialAllCombos, statsTypes = statsTypes, statsNames = statsNames, catLabels = catLabels, b1=b1, b2=b2, mc.cores=mcThreads, mc.preschedule=FALSE))

        calcNCols <- length(unlist(statsDataResQuick[1]))
        statsDataWholeQuick <- matrix(unlist(statsDataResQuick), ncol = calcNCols, byrow=TRUE)
        longCombos <- as.numeric(statsDataWholeQuick[which(statsDataWholeQuick[,7] == -1),6])
        system.time(statsDataResLong <- onlyLongTests(longCombos = longCombos, allCombos = partialAllCombos, statsTypes = statsTypes, statsNames = statsNames, catLabels = catLabels, b1=b1, b2=b2, pThreads=pThreads))
        statsDataResShort <- statsDataWholeQuick[which(statsDataWholeQuick[,7] != -1),]
        statsDataResWhole <- rbind(statsDataResShort,statsDataResLong)
        statsDataWhole <- statsDataResWhole[,-6]

    }else{
        pThreads <- nCores
        mcThreads <- nCores

        outputLine <- paste(Sys.time(), "-- Reached phase 1/2, isolating potentially long running tests", sep=" ")
        print(outputLine)
        ### stderr(outputLine)

        system.time(statsDataResQuick <- mclapply(1:dim(allCombos)[1], onlyQuickTests, allCombos = allCombos, statsTypes = statsTypes, statsNames = statsNames, catLabels = catLabels, tailTypes=tailTypes, b1=b1, b2=b2, writeDir=writeDir, mc.cores=mcThreads, mc.preschedule=FALSE))

        calcNCols <- length(unlist(statsDataResQuick[1]))
        statsDataWholeQuick <- matrix(unlist(statsDataResQuick), ncol = calcNCols, byrow=TRUE)
        longCombos <- as.numeric(statsDataWholeQuick[which(statsDataWholeQuick[,7] == -1),6])
        save.image(paste(writeDir, "/individualRepeatTests/workspaceImage-intermediate-1-", repType, ".RData", sep = ""))

        outputLine <- paste(Sys.time(), "-- Reached phase 2/2, will speed up once past the repeats with high counts", sep=" ")
        print(outputLine)
        ### stderr(outputLine)

        system.time(statsDataResLong <- onlyLongTests(longCombos = longCombos, allCombos = allCombos, statsTypes = statsTypes, statsNames = statsNames, catLabels = catLabels, tailTypes=tailTypes, b1=b1, b2=b2, writeDir=writeDir, pThreads=pThreads))
        save.image(paste(writeDir, "/individualRepeatTests/workspaceImage-intermediate-2-", repType, ".RData", sep = ""))
        statsDataResShort <- statsDataWholeQuick[which(statsDataWholeQuick[,7] != -1),]
        statsDataResWhole <- rbind(statsDataResShort,statsDataResLong)
        statsDataWhole <- statsDataResWhole[,-6]
        #save.image(paste(writeDir, "/individualRepeatTests/workspaceImage-merged-", repType, ".RData", sep = ""))
    }

    save.image(paste(writeDir, "/individualRepeatTests/workspaceImage.RData", sep = ""))

    ### resultColNames <- c("unique_name", "test_type", "repeat_level", "repeat", b1colName, b2colName, pColNames)
    resultColNames <- c("unique_name", "repeat_name", "character_tested", "repeat_classification", b1colName, b2colName)
    numResultColumns <- length(resultColNames) + 2 + sum(tailSelection)
    fullRawResults <- matrix(, ncol=numResultColumns)
    
    for ( i in 1:3){
        repType <- catLabels[i]
        statsData <- statsDataWhole[which(as.integer(statsDataWhole[,2]) == i),]

        ### statsDataCounts <- cbind(paste(repType, statsData[,3], statsTypes[as.integer(statsData[,1])], sep="__"), statsData)

	### statsInfo <- cbind(paste(repType, statsData[,3], statsTypes[as.integer(statsData[,1])], sep="__"), statsData[,1:5]) 
        statsInfo <- cbind(paste(repType, statsData[,3], statsTypes[as.integer(statsData[,1])], sep="__"), statsData[,3], statsTypes[as.integer(statsData[,1])], repType, statsData[,4:5]) 
        colnames(statsInfo) <- resultColNames
        rownames(statsInfo) <- statsInfo[,1]

        statsPData <- statsData[,-c(1,2,3,4,5)]
        statsPData <- statsPData[,as.logical(tailSelection)]
        colnames(statsPData) <- pColNames
        rownames(statsPData) <- statsInfo[,1]

        tenthOfMax <- as.integer(max(max(as.integer(statsData[,4])), max(as.integer(statsData[,4])))/10)
        empCut <- cbind(rep(tenthOfMax, dim(statsInfo)[1]), rep(as.integer(tenthOfMax/10), dim(statsInfo)[1]))
        colnames(empCut) <- c("tenth of max count for repType", "hundredth of max count for repType")
        sectionResults <- cbind(statsInfo, empCut, statsPData)

        fullRawResults <- rbind(fullRawResults, sectionResults)
    }   

    fullRawResults <- fullRawResults[which(!(is.na(fullRawResults[,2]))), ]
    fullRawResultsDF <- as.data.frame(fullRawResults)

    for ( i in 1:3){
        repType <- catLabels[i]
        print(paste(Sys.time(), "-- adjusting the pvalues for grouping:", repType, sep=" "))
        partialRawResultsDF <- fullRawResultsDF[which(fullRawResultsDF$repeat_classification == repType),]
        partialAdjResults <- adjustPValues(repeatDF = partialRawResultsDF, pColNames = pColNames)
        write.csv(lapply(partialAdjResults,as.character), file = paste(writeDir, "/individualRepeatTests/combinedStats--all--", repType, ".csv", sep = ""))
    
        for( st in iTests[iMetrics]) {###this is for the columns in the data files
            print(paste(Sys.time(), "-- adjusting the pvalues for test type:", statsTypes[st], sep=" "))
            partialPartialRawResultsDF <- as.data.frame(partialRawResultsDF[which(partialRawResultsDF$character_tested == statsTypes[st]),])
            partialPartialAdjResults <- adjustPValues(repeatDF = partialPartialRawResultsDF, pColNames = pColNames)
            write.csv(lapply(partialPartialAdjResults,as.character), file = paste(writeDir, "/individualRepeatTests/", repType,"/summaryStats", "--", repType, "--", statsTypes[st], "--", repType, ".csv", sep = ""))
        }
    }
}

print(paste(Sys.time(), "-- adjusting the pvalues for all of the data", sep=" "))
fullAdjResults <- adjustPValues(repeatDF = fullRawResultsDF, pColNames = pColNames)
#write.table(fullAdjResults, file=paste(writeDir, "/individualRepeatTests/combinedStats--all--all.csv", sep=""), col.names=T, row.names=F, sep=",", quote=FALSE)
write.csv(fullAdjResults, file = paste(writeDir, "/individualRepeatTests/combinedStats--all--all.csv", sep = ""))

save.image(paste(writeDir, "/completed_workspace.RData", sep = ""))
resultTmpDir <- paste(writeDir, "/individualRepeatTests", sep="")
resultFinalDir <- paste("./", sep="")

';


if($runR){
    print "analyzing repeats now\n";
    #print "tail -f ". $wd ."/stats.R.out\n";
    #system("R CMD BATCH $rFile ". $rFile .".out");
    system("Rscript $rFile ");
    
    if(!$debug){
	print "cleaning up working directory\n";
	system("mkdir -p ". $wd ."/rawData");
	system("mv ". $wd ."/\*\.tsv ". $wd ."/rawData");
	system("mv ". $wd ."/\*\.bed ". $wd ."/rawData");
	system("mv ". $wd ."/\*\.genome ". $wd ."/rawData");
	system("mv ". $rFile ." ". $wd ."/rawData/");
	system("tar -czf ". $wd ."/rawData.tar.gz ". $wd ."/rawData \&\& rm -rf ". $wd ."/rawData");
    }	
}

print "Completed\n";
