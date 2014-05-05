#!/usr/local/bin/perl

###Open input files

open (FM, "imputed_founder_chrs_2.txt") || die;


###Load input files into memory
##Read in founder matrix data data
print STDOUT "Reading in founder chromosome data.\n";

#Get header for founder matrix and make arrays to hold data
$SEQdataheader = <FM>;
chop $SEQdataheader;
@toks2 = split('\t', $SEQdataheader);
$SEQdataarraylength = scalar (@toks2);
#print STDOUT "Seq data array length is $SEQdataarraylength.\n";
#make array for each column
foreach $colname (@toks2){
	$currentarray{$colname} = [];
	
}

#Read in data and store in the arrays created
$SEQdatalinecount = 0;
while($SEQdataline = <FM>){
	chop $SEQdataline;
	@toks3 = split('\t', $SEQdataline);
	$SEQdatacolcount = 0;
	foreach $column (@toks3){
		$currentarrayarray = @toks2[$SEQdatacolcount];
		push @{$currentarray{$currentarrayarray}}, $column;
		$SEQdatacolcount++;
	}
	$SEQdatalinecount++;
}
print "There are $SEQdatalinecount markers to impute.\n";
$numberofsequencedindividuals = $SEQdataarraylength - 3;
print "There are $numberofsequencedindividuals founding chromosomes.\n";
close (SEQdata);

$founderslistholder = 0;
for ($i=3; $i < $SEQdataarraylength; $i++){
	@toks8 = split('_',@toks2[$i]);
	print "@toks8[0]\n";
	@founderslist1[$founderslistholder] = @toks8[0];
	$founderslistholder++;
}

my @uniquefounders;
my %seenfounders;

foreach my $value (@founderslist1) {
	if (! $seenfounders{$value}++ ) {
		push @uniquefounders, $value;
	}
}

###Create an array of arrays to store the diploid genotypes for each founder
foreach $founderchrforarray (@uniquefounderchrs){
	$arrayoffounderchrs{$founderchrforarray} = [];
}






foreach $currentfounder (@uniquefounders){
	print STDOUT "Constructing diploid genotypes for $currentfounder.\n";
	for($i = 0; $i < $SEQdatalinecount; $i++){
		$sum = ();
		$founder1name = "";
		$founder2name = "";
		$founder1 = ();
		$founder2 = ();
		$founder1name = "$currentfounder"."_1";
		$founder1 = $currentarray{$founder1name}[$i];
		$founder2name = "$currentfounder"."_2";
		$founder2 = $currentarray{$founder2name}[$i];
		if(($founder1 =~ /[0-9]/) && ($founder2 =~ /[0-9]/)) {
			$sum = $founder1 + $founder2;
			$arrayoffounderchrs{$currentfounder}[$i] = $sum;
		}
		else{
			$arrayoffounderchrs{$currentfounder}[$i] = "-";
		}
		#print STDOUT "$founder1name\t$founder1\t$founder2name\t$founder2\tsum\t$sum\n";
	}
}

print STDOUT "Finished imputing diploid genotypes.\nWriting output to GAW_diploid_matrix.txt.\n";
open (OUTFILE, ">GAW_founder_diploid_matrix.txt");

print OUTFILE "@toks2[0]\t@toks2[1]\t@toks2[2]";
foreach $currentfounder (@uniquefounders){
	print OUTFILE "\t$currentfounder";
}
print OUTFILE "\n";

$positionname = @toks2[2];
$snpname = @toks2[0];
$chrname = @toks2[1];

for($i = 0; $i < $SEQdatalinecount; $i++){
	print OUTFILE "$currentarray{$snpname}[$i]\t$currentarray{$chrname}[$i]\t$currentarray{$positionname}[$i]";
	foreach $currentfounder (@uniquefounders){
		print OUTFILE "\t$arrayoffounderchrs{$currentfounder}[$i]";
	}
	print OUTFILE "\n";
}





