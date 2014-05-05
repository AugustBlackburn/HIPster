#!/usr/local/bin/perl

###Open input files

$sequencedata = $ARGV[0];

open (Imputed, "GAW_imputed_diploid_matrix.txt") || die;
open (Sequenced, $sequencedata) || die;
open (OUTFILE, ">GAW_combined_sequenced_and_imputed_matrix.txt") || die;

$Imputedheaderline = <Imputed>;
chop $Imputedheaderline;
@toks1 = split('\t', $Imputedheaderline);
$toks1length = scalar (@toks1);

$Sequencedheaderline = <Sequenced>;
chop $Sequencedheaderline;
@toks2 = split('\t', $Sequencedheaderline);
$toks2length = scalar (@toks2);

foreach $Imputedcolname (@toks1){
	$Imputedarray{$Imputedcolname} = [];	
}

foreach $Sequencedcolname (@toks2){
	$Sequencedarray{$Sequencedcolname} = [];
	
}

my @indivs;
my %seenindivs;
$uniqueindivs = 0;
for ($i=3; $i < $toks1length; $i++) {
	$value = @toks1[$i];
	if (! $seenindivs{$value}++ ) {
		push @indivs, $value;
		$uniqueindivs++;
		print STDOUT "$value\n";
	}
}
for ($i=3; $i < $toks2length; $i++) {
	$value = @toks2[$i];
	if (! $seenindivs{$value}++ ) {
		push @indivs, $value;
		$uniqueindivs++;
		print STDOUT "$value\n";
	}
}

print STDOUT "There are $uniqueindivs people with diploid genotypes.\n";

foreach $diploidindiv(@indivs){
	$diploidarray{$diploidindiv} = [];
}
print OUTFILE "@toks1[0]\t@toks1[1]\t@toks1[2]";

for ($i=0; $i < $uniqueindivs; $i++){
	print OUTFILE "\t@indivs[$i]";
}
print OUTFILE "\n";

while ($line = <Imputed>){
	$line2 = <Sequenced>;
	chop $line;
	chop $line2;
	@toks3 = split('\t', $line);
	@toks4 = split('\t', $line2);
	$Imputedcolcount = 0;
	foreach $column (@toks3){
		$currentimputedarray = @toks1[$Imputedcolcount];
		$Imputedarray{$currentimputedarray}[0]=$column;
		$Imputedcolcount++;
	}
	$Imputedcolcount2 = 0;
	foreach $column2 (@toks4){
		$currentsequencedarray = @toks2[$Imputedcolcount2];
		$Sequencedarray{$currentsequencedarray}[0]=$column2;
		$Imputedcolcount2++;
	}
	foreach $diploiddude(@indivs){
		$Genotype1="";
		$Genotype2="";
		$Decidedgenotype=();
		$Genotype1=$Imputedarray{$diploiddude}[0];
		$Genotype2=$Sequencedarray{$diploiddude}[0];
		if($Genotype2 =~ /[0-9]/){
			$Decidedgenotype=$Genotype2;
		}
		elsif($Genotype1 =~ /[0-9]/){
			$Decidedgenotype=$Genotype1;
		}
		else{
			$Decidedgenotype="-";
		}
		$diploidarray{$diploiddude}[0] = $Decidedgenotype;
	}
print OUTFILE "@toks3[0]\t@toks3[1]\t@toks3[2]";
for ($i=0; $i < $uniqueindivs; $i++){
	$uniqueid = @indivs[$i];
	print OUTFILE "\t$diploidarray{$uniqueid}[0]";
}
print OUTFILE "\n";
}
