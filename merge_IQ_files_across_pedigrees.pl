#!/usr/local/bin/perl

###Open input files

open (IQS2, "IQS_2.txt") || die;
open (IQS3, "IQS_3.txt") || die;
open (IQS4, "IQS_4.txt") || die;
open (IQS5, "IQS_5.txt") || die;
open (IQS6, "IQS_6.txt") || die;
open (IQS7, "IQS_7.txt") || die;
open (IQS8, "IQS_8.txt") || die;
open (IQS9, "IQS_9.txt") || die;
open (IQS10, "IQS_10.txt") || die;
open (IQS11, "IQS_11.txt") || die;
open (IQS14, "IQS_14.txt") || die;
open (IQS15, "IQS_15.txt") || die;
open (IQS16, "IQS_16.txt") || die;
open (IQS17, "IQS_17.txt") || die;
open (IQS20, "IQS_20.txt") || die;
open (IQS21, "IQS_21.txt") || die;
open (IQS23, "IQS_23.txt") || die;
open (IQS25, "IQS_25.txt") || die;
open (IQS27, "IQS_27.txt") || die;
open (IQS47, "IQS_47.txt") || die;

open (OUTFILE, ">IQS_all.txt") || die;

$Headerline2 = <IQS2>;
$Headerline3 = <IQS3>;
$Headerline4 = <IQS4>;
$Headerline5 = <IQS5>;
$Headerline6 = <IQS6>;
$Headerline7 = <IQS7>;
$Headerline8 = <IQS8>;
$Headerline9 = <IQS9>;
$Headerline10 = <IQS10>;
$Headerline11 = <IQS11>;
$Headerline14 = <IQS14>;
$Headerline15 = <IQS15>;
$Headerline16 = <IQS16>;
$Headerline17 = <IQS17>;
$Headerline20 = <IQS20>;
$Headerline21 = <IQS21>;
$Headerline23 = <IQS23>;
$Headerline25 = <IQS25>;
$Headerline27 = <IQS27>;
$Headerline47 = <IQS47>;
chop $Headerline2;

print OUTFILE "$Headerline2\n";

while ($line2 = <IQS2>) {
	$line3 = <IQS3>;
	$line4 = <IQS4>;
	$line5 = <IQS5>;
	$line6 = <IQS6>;
	$line7 = <IQS7>;
	$line8 = <IQS8>;
	$line9 = <IQS9>;
	$line10 = <IQS10>;
	$line11 = <IQS11>;
	$line14 = <IQS14>;
	$line15 = <IQS15>;
	$line16 = <IQS16>;
	$line17 = <IQS17>;
	$line20 = <IQS20>;
	$line21 = <IQS21>;
	$line23 = <IQS23>;
	$line25 = <IQS25>;
	$line27 = <IQS27>;
	$line47 = <IQS47>;
	chop $line2;
	chop $line3;
	chop $line4;
	chop $line5;
	chop $line6;
	chop $line7;
	chop $line8;
	chop $line9;
	chop $line10;
	chop $line11;
	chop $line14;
	chop $line15;
	chop $line16;
	chop $line17;
	chop $line20;
	chop $line21;
	chop $line23;
	chop $line25;
	chop $line27;
	chop $line47;
	@toks2 = split('\t', $line2);
	@toks3 = split('\t', $line3);
	@toks4 = split('\t', $line4);
	@toks5 = split('\t', $line5);
	@toks6 = split('\t', $line6);
	@toks7 = split('\t', $line7);
	@toks8 = split('\t', $line8);
	@toks9 = split('\t', $line9);
	@toks10 = split('\t', $line10);
	@toks11 = split('\t', $line11);
	@toks14 = split('\t', $line14);
	@toks15 = split('\t', $line15);
	@toks16 = split('\t', $line16);
	@toks17 = split('\t', $line17);
	@toks20 = split('\t', $line20);
	@toks21 = split('\t', $line21);
	@toks23 = split('\t', $line23);
	@toks25 = split('\t', $line25);
	@toks27 = split('\t', $line27);
	@toks47 = split('\t', $line47);
	$Totaln = @toks2[3] + @toks3[3] + @toks4[3] + @toks5[3] + @toks6[3] + @toks7[3] + @toks8[3] + @toks9[3] + @toks10[3] + @toks11[3] + @toks14[3] + @toks15[3] + @toks16[3] + @toks17[3] + @toks20[3] + @toks21[3] + @toks23[3] + @toks25[3] + @toks27[3] + @toks47[3];
	$Total0 = @toks2[4] + @toks3[4] + @toks4[4] + @toks5[4] + @toks6[4] + @toks7[4] + @toks8[4] + @toks9[4] + @toks10[4] + @toks11[4] + @toks14[4] + @toks15[4] + @toks16[4] + @toks17[4] + @toks20[4] + @toks21[4] + @toks23[4] + @toks25[4] + @toks27[4] + @toks47[4];
	$Total1 = @toks2[5] + @toks3[5] + @toks4[5] + @toks5[5] + @toks6[5] + @toks7[5] + @toks8[5] + @toks9[5] + @toks10[5] + @toks11[5] + @toks14[5] + @toks15[5] + @toks16[5] + @toks17[5] + @toks20[5] + @toks21[5] + @toks23[5] + @toks25[5] + @toks27[5] + @toks47[5];
	$Total2 = @toks2[6] + @toks3[6] + @toks4[6] + @toks5[6] + @toks6[6] + @toks7[6] + @toks8[6] + @toks9[6] + @toks10[6] + @toks11[6] + @toks14[6] + @toks15[6] + @toks16[6] + @toks17[6] + @toks20[6] + @toks21[6] + @toks23[6] + @toks25[6] + @toks27[6] + @toks47[6];
	$Totalaccurate = @toks2[7] + @toks3[7] + @toks4[7] + @toks5[7] + @toks6[7] + @toks7[7] + @toks8[7] + @toks9[7] + @toks10[7] + @toks11[7] + @toks14[7] + @toks15[7] + @toks16[7] + @toks17[7] + @toks20[7] + @toks21[7] + @toks23[7] + @toks25[7] + @toks27[7] + @toks47[7];
	if($Totaln == 0 ){
		$Concordance = "NA";
		$ByChance = "NA";
	}else{
		$Concordance = $Totalaccurate / $Totaln ;
		$ByChance = (($Total0 * $Total0) + ($Total1 * $Total1) + ($Total2 * $Total2))/($Totaln * $Totaln);
	}
	if(($ByChance || $Concordance) eq "NA"){
		$IQS = "NA";
	}elsif((1 - $ByChance) == 0){
		$IQS = "NA";
	}else{
		$IQS = ($Concordance - $ByChance)/(1 - $ByChance);
	}
print OUTFILE "@toks2[0]\t@toks2[1]\t@toks2[2]\t$Totaln\t$Total0\t$Total1\t$Total2\t$Totalaccurate\t$Concordance\t$ByChance\t$IQS\n";
}

close OUTFILE;
close IQS2;
close IQS3;
close IQS4;
close IQS5;
close IQS6;
close IQS7;
close IQS8;
close IQS9;
close IQS10;
close IQS11;
close IQS14;
close IQS15;
close IQS16;
close IQS17;
close IQS20;
close IQS21;
close IQS23;
close IQS25;
close IQS27;
close IQS47;
