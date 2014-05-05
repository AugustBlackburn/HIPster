#!/usr/local/bin/perl

#$ARGV[0] is the converted 012 file
open (INFILE_GENO, $ARGV[0]) || die;
open (OUTFILE1, ">First_half_GWAS.txt") || die;
open (OUTFILE2, ">Second_half_GWAS.txt") || die;

$headerline = <INFILE_GENO>;
chop $headerline;
print OUTFILE1 "$headerline\n";
print OUTFILE2 "$headerline\n";


$counter = 1;
while($line = <INFILE_GENO>){
	if($counter == 1){
		chop $line;
		print OUTFILE1 "$line\n";
		$counter = 0;
	}elsif($counter == 0){
		chop $line;
		print OUTFILE2 "$line\n";
		$counter = 1;
	}else{
		print STDOUT "Error.\n";
	}
}
close INFILE_GENO;
close OUTFILE1;
close OUTFILE2;

