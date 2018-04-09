#!/usr/bin/perl

if(scalar(@ARGV)<1) {
	print "perl cv_partition.pl 0\n";
	die;
}

$path0="../../../data/trimmed_set/";
$path1="../../../data/normalization/raw/";
$path2="../../../data/cv_set/";
$path3="./";

%rna=();
open RNASEQ, "${path1}trimmed_breast_rna.txt.avg_imputation" or die; # read rna-seq
$line=<RNASEQ>;
chomp $line;
@ids=split "\t", $line;
shift @ids;
while($line=<RNASEQ>){
	chomp $line;
	@vals=split "\t", $line;
	$gene=shift @vals;
	$i=0;
	foreach $sample (@ids){
		$rna{$sample}{$gene}=$vals[$i];
		$i++;
	}
}
close RNASEQ;

%cnv=();
open CNV, "${path1}trimmed_breast_cnv.txt.avg_imputation" or die; # read rna-seq
$line=<CNV>;
chomp $line;
@ids=split "\t", $line;
shift @ids;
while($line=<CNV>){
	chomp $line;
	@vals=split "\t", $line;
	$gene=shift @vals;
	$i=0;
	foreach $sample (@ids){
		$cnv{$sample}{$gene}=$vals[$i];
		$i++;
	}
}
close CNV;

foreach $num (@ARGV){
	open TEST_P, "${path2}test_breast_proteome_${num}.txt" or die;	
	$line=<TEST_P>;
	chomp $line;
	@test_ids=split "\t", $line;
	shift @test_ids;

	open OUTPUT_PRED, ">${path3}prediction_breast_proteome_${num}.txt" or die;
	print OUTPUT_PRED "$line\n";

	while($line=<TEST_P>){
		chomp $line;
		@vals=split "\t", $line;
		$gene=shift @vals;

		print OUTPUT_PRED "$gene";
		foreach $id (@test_ids){
			$tmp=$rna{$id}{$gene};
			print OUTPUT_PRED "\t$tmp";
		}
		print OUTPUT_PRED "\n";
	}
	close TEST_P;
	close OUTPUT_PRED;
}

