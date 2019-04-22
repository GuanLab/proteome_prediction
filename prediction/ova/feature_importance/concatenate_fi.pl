#!/usr/bin/perl

## Here I concatenate the separated feature importance file into one

($list_gene) = @ARGV;

open OUTPUT, ">fi.txt" or die;
open INPUT, "$list_gene" or die; # header
print OUTPUT "FI";
while($line=<INPUT>){
    chomp $line;
    print OUTPUT "\t$line";
}
print OUTPUT "\n";
close INPUT;


@list_file=glob ("importance/*txt");
foreach $file (@list_file){
    @tmp=split '/', $file;
    @tmp1=split '\.', $tmp[1]; # You need to use single quote here! Be careful!
    $gene=$tmp1[0];
    print OUTPUT "$gene";
    open INPUT, "${file}" or die;
    while($line=<INPUT>){
        chomp $line;
        print OUTPUT "\t$line";
    }
    print OUTPUT "\n";
    close INPUT;
}

