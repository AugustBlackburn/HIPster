#!/usr/local/bin/perl

###Open input files

open (FM, "GAW_founder_diploid_matrix.txt") || die;
open (NFM, "GAW_diploid_matrix.txt") || die;
open (OUTFILE, ">GAW_imputed_diploid_matrix.txt") || die;

$FMheaderline = <FM>;
chop $FMheaderline;
@toks1 = split('\t', $FMheaderline);
$toks1length = scalar (@toks1);

$NFMheaderline = <NFM>;
chop $NFMheaderline;
@toks2 = split('\t', $NFMheaderline);
$toks2length = scalar (@toks2);

print OUTFILE "@toks1[0]";
for($i = 1; $i < $toks1length; $i++){
	print OUTFILE "\t@toks1[$i]";
}
for($j = 3; $j < $toks2length; $j++){
	print OUTFILE "\t@toks2[$j]";
}
print OUTFILE "\n";

while ($line = <FM>){
	$line2 = <NFM>;
	chop $line;
	chop $line2;
	@toks3 = split('\t', $line);
	@toks4 = split('\t', $line2);
	print OUTFILE "@toks3[0]";
	for($i = 1; $i < $toks1length; $i++){
		print OUTFILE "\t@toks3[$i]";
	}
	for($j = 3; $j < $toks2length; $j++){
		print OUTFILE "\t@toks4[$j]";
	}
	print OUTFILE "\n";
}
