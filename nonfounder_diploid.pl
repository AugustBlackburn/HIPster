#!/usr/local/bin/perl

###Open input files

open (FM, "Imputed_founder_chrs.txt") || die;
open (NFM, "GAW_nonfounder_matrix.txt") || die;

###Load input files into memory
##Read in non-founder matrix
print STDOUT "Reading in non-founder matrix.\n";
$NFMheader = <NFM>;
$NFMlinecount = 0;
while($NFMline = <NFM>){
	chop $NFMline;
	@toks = split('\t', $NFMline);
	@NFMindiv[$NFMlinecount] = $toks[0];
	@NFMfoundingchr[$NFMlinecount] = $toks[1];
	@NFMstart[$NFMlinecount] = $toks[2];
	@NFMend[$NFMlinecount] = $toks[3];
	@NFMPath[$NFMlinecount] = $toks [4];
	$NFMlinecount++;
}
close (NFM);

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

###This is code for checking to make sure the data was read in correctly.
$NFMarraylength = scalar(@NFMindiv);
#print STDOUT "NFM array length is $NFMarraylength.\n";
#for($i = 0; $i < $NFMarraylength; $i++){
#	print STDOUT "@NFMindiv[$i],@NFMfoundingchr[$i],@NFMstart[$i],@NFMend[$i],@NFMPath[$i]\n";
#}


#for( $i = 0; $i < 50; $i++) {
#	foreach $colname (@toks2){
#		print "$currentarray{$colname}[$i]\t";
#	}
#	print "\n";
#}



###Identify all of the non-founders
my @nonfounders;
my %seennonfounders;
$uniqueNFindivids = 0;
print STDOUT "The non-founders are:\n";
foreach my $value (@NFMindiv) {
	if (! $seennonfounders{$value}++ ) {
		push @nonfounders, $value;
		$uniqueNFindivids++;
		print STDOUT "$value\n";
	}
}
print STDOUT "There are $uniqueNFindivids non-founders.\n";

###Create an array of arrays to store the diploid genotypes for each non-founder
foreach $nonfounderforarray (@nonfounders){
	$arrayofnonfounders{$nonfounderforarray} = [];
}

###Universal variables
$positionname = @toks2[2];
$snpname = @toks2[0];
$chrname = @toks2[1];

foreach $currentnonfounder (@nonfounders){
	print STDOUT "Constructing diploid genotypes for $currentnonfounder.\n";
	@currentNFMindiv = [];
	@currentNFMfoundingchr = [];
	@currentNFMstart = [];
	@currentNFMend = [];
	@currentNFMPath = [];	
	$currentNFMlinecount = 0;
	for($i = 0; $i < $NFMarraylength; $i++) {
		if (@NFMindiv[$i] eq $currentnonfounder){
			#print STDOUT "@NFMfoundingchr[$i] = $currentfoundingchromosome\n";
			@currentNFMindiv[$currentNFMlinecount] = @NFMindiv[$i];
			@currentNFMfoundingchr[$currentNFMlinecount] = @NFMfoundingchr[$i];
			@currentNFMstart[$currentNFMlinecount] = @NFMstart[$i];
			@currentNFMend[$currentNFMlinecount] = @NFMend[$i];
			@currentNFMPath[$currentNFMlinecount] = @NFMPath[$i];
			$currentNFMlinecount++;
		}
	}
	$currentNFMarraylength = scalar(@currentNFMindiv);
	#print STDOUT "The current NFMarray length is $currentNFMarraylength.\n";
	#for($i = 0; $i < $currentNFMarraylength; $i++) {
	#	print "@currentNFMindiv[$i]\t@currentNFMfoundingchr[$i]\t@currentNFMstart[$i]\t@currentNFMend[$i]\t@currentNFMPath[$i]\n";
	#}
	#print "\n";
	
	for($i = 0; $i < $SEQdatalinecount; $i++){
		$sum = ();
		$founder1 = "";
		$founder2 = "";
		$has_f = 0;
		
		for($j = 0; $j < $currentNFMarraylength; $j++){
			if (($currentarray{$positionname}[$i] >= @currentNFMstart[$j]) && ($currentarray{$positionname}[$i] <= @currentNFMend[$j])){
				if($has_f == 0){
					$founder1 = @currentNFMfoundingchr[$j];
					$has_f = 1;
				}
				else{
					$founder2 = @currentNFMfoundingchr[$j];
					$has_f = 2;
				}
			}
		}
		if ($has_f <2){
			$arrayofnonfounders{$currentnonfounder}[$i] = "-";
		}
		else{
			if(($currentarray{$founder1}[$i] =~ /[0-9]/) && ($currentarray{$founder2}[$i] =~ /[0-9]/)) {
				$sum = $currentarray{$founder1}[$i] + $currentarray{$founder2}[$i];
				$arrayofnonfounders{$currentnonfounder}[$i] = $sum;
			}
			else{
				$arrayofnonfounders{$currentnonfounder}[$i] = "-";
			}
		}

	}
}

print STDOUT "Finished imputing diploid genotypes.\nWriting output to GAW_diploid_matrix.txt.\n";
open (OUTFILE, ">GAW_diploid_matrix.txt");

print OUTFILE "@toks2[0]\t@toks2[1]\t@toks2[2]";
foreach $currentnonfounder (@nonfounders){
	print OUTFILE "\t$currentnonfounder";
}
print OUTFILE "\n";
for($i = 0; $i < $SEQdatalinecount; $i++){
	print OUTFILE "$currentarray{$snpname}[$i]\t$currentarray{$chrname}[$i]\t$currentarray{$positionname}[$i]";
	foreach $currentnonfounder (@nonfounders){
		print OUTFILE "\t$arrayofnonfounders{$currentnonfounder}[$i]";
	}
	print OUTFILE "\n";
}
