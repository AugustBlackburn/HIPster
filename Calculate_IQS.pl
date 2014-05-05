#!/usr/local/bin/perl

###Open input files

$sequencedata = $ARGV[0];

open (Imputed, "GAW_combined_sequenced_and_imputed_matrix.txt") || die;
open (Sequenced, $sequencedata) || die;
open (OUTFILE, ">Checking_imputation_accuracy_by_marker.txt") || die;
open (ID_OUTFILE, ">Checking_imputation_accuracy_by_id.txt") || die;
open (OUTFILE3_IQS, ">IQS_statistics.txt") || die;

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
print OUTFILE "@toks1[0]\t@toks1[1]\t@toks1[2]\tAccurately_typed\tAccurately_typed_hetero\tInaccurately_typed\tInaccurately_typed_hetero\tImputation_added_homozygous\tImputation_added_heterozygous\tNot_imputed_homozygous\tNot_imputed_heterozygous\n";



foreach $diploiddude(@indivs){
	$accuratecount{$diploiddude}= 0;
	$inaccuratecount{$diploiddude}= 0;
	$accurate1count{$diploiddude}= 0;
	$inaccurate1count{$diploiddude}= 0;
	$notimputedhomo{$diploiddude}= 0;
	$imputedhomo{$diploiddude}= 0;
	$imputedhetero{$diploiddude}= 0;
	$notimputedhetero{$diploiddude}= 0;		
}

$accuratecount{$marker}= 0;
$inaccuratecount{$marker}= 0;
$accurate1count{$marker}= 0;
$inaccurate1count{$marker}= 0;
$notimputedhomo{$marker}= 0;
$imputedhomo{$marker}= 0;
$imputedhetero{$marker}= 0;
$notimputedhetero{$marker}= 0;		


while ($line = <Imputed>){
	$line2 = <Sequenced>;
	chop $line;
	chop $line2;
	@toks3 = split('\t', $line);
	@toks4 = split('\t', $line2);
	$Imputedcolcount = 0;
	
	$marker = @toks3[0];
	$accuratecount{$marker}= 0;
	$inaccuratecount{$marker}= 0;
	$accurate1count{$marker}= 0;
	$inaccurate1count{$marker}= 0;
	$notimputedhomo{$marker}= 0;
	$imputedhomo{$marker}= 0;
	$imputedhetero{$marker}= 0;
	$notimputedhetero{$marker}= 0;
	push(@markers,$marker);
	$local1{$marker}=@toks3[1];
	$local2{$marker}=@toks3[2];
	$accurate00{$marker} = 0;
	$accurate11{$marker} = 0;
	$accurate22{$marker} = 0;
	$inaccuratecount01{$marker} = 0;
	$inaccuratecount02{$marker} = 0;
	$inaccuratecount10{$marker} = 0;
	$inaccuratecount12{$marker} = 0;
	$inaccuratecount20{$marker} = 0;
	$inaccuratecount21{$marker} = 0;
	
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
		if(($Genotype2 =~ /[0-9]/)&&($Genotype1 =~ /[0-9]/)){
			if(($Genotype2 =~ /0/)&&($Genotype1 =~ /0/)){
				$accuratecount{$diploiddude}++;
				$accuratecount{$marker}++;
				$accurate00{$marker}++;
			}
			elsif(($Genotype2 =~ /1/)&&($Genotype1 =~ /1/)){
				$accuratecount{$diploiddude}++;
				$accurate1count{$diploiddude}++;	
				$accuratecount{$marker}++;
				$accurate1count{$marker}++;
				$accurate11{$marker}++;				
			
			}
			elsif(($Genotype2 =~ /2/)&&($Genotype1 =~ /2/)){
				$accuratecount{$diploiddude}++;
				$accuratecount{$marker}++;
				$accurate22{$marker}++;
			}
			elsif(($Genotype2 =~ /0/)&&($Genotype1 =~ /1/)){
				$inaccuratecount{$diploiddude}++;
				$inaccuratecount{$marker}++;
				$inaccuratecount01{$marker}++;
			}
			elsif(($Genotype2 =~ /0/)&&($Genotype1 =~ /2/)){
				$inaccuratecount{$diploiddude}++;
				$inaccuratecount{$marker}++;
				$inaccuratecount02{$marker}++;
			}
			elsif(($Genotype2 =~ /1/)&&($Genotype1 =~ /0/)){
				$inaccurate1count{$diploiddude}++;
				$inaccuratecount{$diploiddude}++;
				$inaccurate1count{$marker}++;
				$inaccuratecount{$marker}++;
				$inaccuratecount10{$marker}++;
			}
			elsif(($Genotype2 =~ /1/)&&($Genotype1 =~ /2/)){
				$inaccurate1count{$diploiddude}++;
				$inaccuratecount{$diploiddude}++;
				$inaccurate1count{$marker}++;
				$inaccuratecount{$marker}++;
				$inaccuratecount12{$marker}++;
			}
			elsif(($Genotype2 =~ /2/)&&($Genotype1 =~ /0/)){
				$inaccuratecount{$diploiddude}++;
				$inaccuratecount{$marker}++;
				$inaccuratecount20{$marker}++;
			}
			elsif(($Genotype2 =~ /2/)&&($Genotype1 =~ /1/)){
				$inaccuratecount{$diploiddude}++;
				$inaccuratecount{$marker}++;
				$inaccuratecount21{$marker}++;
			}
		}
		elsif($Genotype1 =~ /[0-9]/){
			if($Genotype1 =~ /0/){
				$imputedhomo{$diploiddude}++;
				$imputedhomo{$marker}++;
			}
			elsif($Genotype1 =~ /1/){
				$imputedhetero{$diploiddude}++;
				$imputedhetero{$marker}++;
			}
			elsif($Genotype1 =~ /2/){
				$imputedhomo{$diploiddude}++;
				$imputedhomo{$marker}++;
			}
		}
		elsif($Genotype2 =~ /[0-9]/){
			if($Genotype2 =~ /0/){
				$notimputedhomo{$diploiddude}++;
				$notimputedhomo{$marker}++;
			}
			elsif($Genotype2 =~ /1/){
				$notimputedhetero{$diploiddude}++;
				$notimputedhetero{$marker}++;
			}
			elsif($Genotype2 =~ /2/){
				$notimputedhomo{$diploiddude}++;
				$notimputedhomo{$marker}++;
			}
			
		}
		####
	}
}

print ID_OUTFILE "ID\tAccurately_typed\tAccurately_typed_hetero\tInaccurately_typed\tInaccurately_typed_hetero\tImputation_added_homozygous\tImputation_added_heterozygous\tNot_imputed_homozygous\tNot_imputed_heterozygous\n";
foreach $diploiddude(@indivs){
print ID_OUTFILE "$diploiddude\t$accuratecount{$diploiddude}\t$accurate1count{$diploiddude}\t$inaccuratecount{$diploiddude}\t$inaccurate1count{$diploiddude}\t$imputedhomo{$diploiddude}\t$imputedhetero{$diploiddude}\t$notimputedhomo{$diploiddude}\t$notimputedhetero{$diploiddude}\n";
}


foreach $marker(@markers){
	print OUTFILE "$marker\t$local1{$marker}\t$local2{$marker}\t$accuratecount{$marker}\t$accurate1count{$marker}\t$inaccuratecount{$marker}\t$inaccurate1count{$marker}\t$imputedhomo{$marker}\t$imputedhetero{$marker}\t$notimputedhomo{$marker}\t$notimputedhetero{$marker}\n";
}
close OUTFILE;
close ID_OUTFILE;
close Imputed;
close Sequenced;

print OUTFILE3_IQS "@toks1[0]\t@toks1[1]\t@toks1[2]\tTotal_n\tTotal_0\tTotal_1\tTotal_2\tTotal_accurate\tConcordance\tBy_chance\tIQS\n";
foreach $marker(@markers){
	
	$totaln = $accurate00{$marker} + $accurate11{$marker} + $accurate22{$marker} + $inaccuratecount01{$marker} + $inaccuratecount02{$marker} + $inaccuratecount10{$marker} + $inaccuratecount12{$marker} + $inaccuratecount20{$marker} + $inaccuratecount21{$marker};
	$total0 = $accurate00{$marker} + $inaccuratecount01{$marker} + $inaccuratecount02{$marker};
	$total1 = $accurate11{$marker} + $inaccuratecount10{$marker} + $inaccuratecount12{$marker};
	$total2 = $accurate22{$marker} + $inaccuratecount20{$marker} + $inaccuratecount21{$marker};
	$totalaccurate = $accurate00{$marker} + $accurate11{$marker} + $accurate22{$marker};
	if(($totaln * $totaln) == 0){
		$bychance = "NA";
	}else{
		$bychance = (($total0 * $total0) + ($total1 * $total1) + ($total2 * $total2))/($totaln * $totaln);
	}
	if($totaln == 0 ){
		$concordance = "NA";
	}else{
		$concordance = $totalaccurate / $totaln ;
	}
	if(($bychance || $concordance) eq "NA"){
		$IQS = "NA";
	}elsif((1 - $bychance) == 0){
		$IQS = "NA";
	}else{
		$IQS = ($concordance - $bychance)/(1 - $bychance);
	}
	print OUTFILE3_IQS "$marker\t$local1{$marker}\t$local2{$marker}\t$totaln\t$total0\t$total1\t$total2\t$totalaccurate\t$concordance\t$bychance\t$IQS\n";
}

	





