#!/usr/local/bin/perl

#Arguments are the matrix of diploid genotypes and the percentage of markers to randomly mask
$SEQMATRIX = $ARGV[0];
$drop = $ARGV[1];

open (SEQdata, $SEQMATRIX) || die;
open (OUTFILE1, ">Masked_genotypes.txt") || die;
open (OUTFILE2, ">Genotypes_not_masked.txt") || die;

##Read in seq data
print STDOUT "Reading in SEQ data.\n";

#Get header for seq data and make arrays to hold data
$SEQdataheader = <SEQdata>;
chop $SEQdataheader;

print OUTFILE1 "$SEQdataheader\n";
print OUTFILE2 "$SEQdataheader\n";
$maskcount = 0;
$keptcount = 0;
while ($line = <SEQdata>){
	chop $line;
	@toks1 = split('\t', $line);
	$toksarraylength = scalar(@toks1);
	print OUTFILE1 "@toks1[0]\t@toks1[1]\t@toks1[2]";
	print OUTFILE2 "@toks1[0]\t@toks1[1]\t@toks1[2]";
	for ($i=3; $i < $toksarraylength; $i++) {
		$value = rand(1);
		if ($value <= $drop){
			print OUTFILE1 "\t@toks1[$i]";
			print OUTFILE2 "\t-";
			$maskcount++;
		}else{
			print OUTFILE2 "\t@toks1[$i]";
			print OUTFILE1 "\t-";
			$keptcount++;
		}
	}
	print OUTFILE1 "\n";
	print OUTFILE2 "\n";
}

close OUTFILE1;
close OUTFILE2;

$maskedpercentage = $maskcount*100/($maskcount+$keptcount);

print STDOUT "$maskedpercentage percentage of genotypes in $SEQMATRIX were masked and stored in Masked_genotypes.txt. The remaining genotypes were written to Genotypes_not_masked.txt\n";

		