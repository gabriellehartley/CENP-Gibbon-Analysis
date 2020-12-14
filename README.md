# Gibbon-RepeatMasker-Script
A script used to summarize RepeatMasker annotations from the .out file.

usage: 

perl ./RepeatMaskerSummarizer.perl file.fa.out

# Repstat 

1. Make bed files from fasta files

	perl ./repstat/makeBed.pl Lib1.fasta
	perl ./repstat/makeBed.pl Lib2.fasta

2. Run RepeatMasker on libraries (independent of Repstat scripts)

3. Combine RepeatMasker .out files of interest

	cat Lib1.fasta.out Lib2.fasta.out > Lib1-2.fasta.out

4. Run Repstat

perl ./repstat-individualTests.pl --rmOut Lib1-2fasta.out --bed1 Lib1.fasta.bed --bed2 Lib2.fasta.bed --outDir Lib1_vs_Lib2_2sided --cores 24

