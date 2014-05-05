#!/usr/local/bin/perl

###Open input files

$sequencedata = $ARGV[0];
$unrelfile = $ARGV[1];

open (Unrel, $unrelfile) || die;
open (Sequenced, $sequencedata) || die;
open (OUTFILE, ">Marker_MAF_calculated_based_on_unrelateds.txt") || die;

$Sequencedheaderline = <Sequenced>;
chop $Sequencedheaderline;
@toks2 = split('\t', $Sequencedheaderline);
$toks2length = scalar (@toks2);

foreach $Sequencedcolname (@toks2){
	$Sequencedarray{$Sequencedcolname} = [];
	
}




my @indivs2;
my %seenindivs2;
$uniqueindivs2 = 0;
while ($line2 = <Unrel>){
	chop $line2;
	@toks4 = split('\t',$line2);
	$value2 = @toks4[0];
	if (! $seenindivs2{$value2}++ ) {
		push @indivs2, $value2;
		$uniqueindivs2++;
		print STDOUT "$value2\n";
	}

}
print STDOUT "MAF to be calculated based on $uniqueindivs2 individuals.\n";



my @indivs;
my %seenindivs;
$uniqueindivs = 0;

for ($i=3; $i < $toks2length; $i++) {
	$value = @toks2[$i];
	if (! $seenindivs{$value}++ ) {
		push @indivs, $value;
		$uniqueindivs++;
		print STDOUT "$value\n";
	}
}

print OUTFILE "@toks2[0]\t@toks2[1]\t@toks2[2]\tMAF\n";

while ($line = <Sequenced>){
	chop $line;
	@toks3 = split('\t', $line);
	$count0 = 0;
	$count1 = 0;
	$count2 = 0;
	$MAF = "NA";	

	$Imputedcolcount2 = 0;
	foreach $column2 (@toks3){
		$currentsequencedarray = @toks2[$Imputedcolcount2];
		$Sequencedarray{$currentsequencedarray}[0]= $column2;
		$Imputedcolcount2++;
	}

	foreach $diploiddude(@indivs2){
		$Genotype2="";
		$Genotype2=$Sequencedarray{$diploiddude}[0];
		if($Genotype2 =~ /[0-9]/){
			if($Genotype2 =~ /[0]/){
				$count0++;
			}elsif($Genotype2 =~ /[1]/){
				$count1++;
			}elsif($Genotype2 =~ /[2]/){
				$count2++;
			}
		}
	}
	if(($count0+$count1+$count2) > 0){
		$MAF = ((2*$count0)+$count1)/(2*($count0+$count1+$count2));	
	}
	if($MAF>0.5){
		$MAF2 = 1-$MAF;
		$MAF = $MAF2;
	}
print OUTFILE "@toks3[0]\t@toks3[1]\t@toks3[2]\t$MAF\n";
}
