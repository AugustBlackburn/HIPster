#!/usr/local/bin/perl

###Open input files

$known = $ARGV[0];
$masked = $ARGV[1];

open (Masked, $masked) || die;
open (Known, $known) || die;
open (OUTFILE, ">Marker_missingness_percentage.txt") || die;

$Knownheaderline = <Known>;
$Maskedheaderline = <Masked>;
chop $Maskedheaderline;
@toks2 = split('\t', $Maskedheaderline);
$toks2length = scalar (@toks2);

print OUTFILE "@toks2[0]\t@toks2[1]\t@toks2[2]\tKnown\tMasked\tPercent_known\n";

while ($line = <Known>){
	$line2 = <Masked>;
	chop $line;
	chop $line2;
	@toks3 = split('\t', $line);
	@toks4 = split('\t', $line2);
	$countKnown = 0;
	$countMasked = 0;
	$PercentKnown = "NA";	

	
	for ($i=3;$i<$toks2length;$i++){
		$Genotype = @toks3[$i];
		$Genotype2 = @toks4[$i];
		if($Genotype =~ /[0-9]/){
			$countKnown++;
		}
		if($Genotype2 =~ /[0-9]/){
			$countMasked++;
		}
	}
	if(($countKnown + $countMasked) > 0){
		$PercentKnown = $countKnown/($countKnown + $countMasked)
 	}
print OUTFILE "@toks3[0]\t@toks3[1]\t@toks3[2]\t$countKnown\t$countMasked\t$PercentKnown\n";
}
