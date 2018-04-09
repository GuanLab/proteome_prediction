#!/usr/bin/perl
#


@all_data=@ARGV;

### take the anchoring point from data set 1
#

open FILE, "$all_data[0]" or die;
<FILE>;
@all=();
while ($line=<FILE>){
	chomp $line;
	@table=split "\t", $line;
	$i=1;
	shift @table;
	foreach $aaa (@table){
		$all[$i].="\t$aaa";
		$i++;
	}
}

close FILE;
	

$imax=$i;
$i=1;
while ($i<$imax){
	@value=split "\t", $all[$i];
	shift @value;
	@value=sort{$a<=>$b}@value;
	$j=0;
	foreach $vvv (@value){
		$sum[$j]+=$vvv;
		$count[$j]++;
		$j++;
	}
	$i++;
}

$jmax=$j;
$j=0;
while ($j<$jmax){
	$avg[$j]=$sum[$j]/$count[$j];
	$j++;
}


foreach $eachdata (@all_data){
	@all=();
	%map=();
	open OLD, "$eachdata" or die;
	$line=<OLD>;
	while ($line=<OLD>){
		chomp $line;
		@table=split "\t", $line;
		$i=1;
		shift @table;
		foreach $aaa (@table){
			$all[$i].="\t$aaa";
			$i++;
		}
	}
	close OLD;
	close NEW;
	
	$imax=$i;
	$i=1;
	while ($i<$imax){
		@value=split "\t", $all[$i];
		shift @value;
		@value=sort{$a<=>$b}@value;
		$jmax_new=scalar(@value);
		$j=0;
		foreach $aaa (@value){
			$new_id=int(($j/$jmax_new)*$jmax);
			$map[$i]{$aaa}=$avg[$new_id];
			$j++;
		}
		$i++;
	}
	open OLD, "$eachdata" or die;
	$new=$eachdata.'.anchor';
	open NEW, ">$new" or die;

	$line=<OLD>;
	print NEW "$line";
	while ($line=<OLD>){
		chomp $line;
		@table=split "\t", $line;
		$i=1;
		$gene=shift @table;
		print NEW "$gene";
		$i=1;
		foreach $aaa (@table){
			print NEW "\t$map[$i]{$aaa}";
			$i++;
		}
		print NEW "\n";
	}
	close OLD;
	close NEW;
}

		



