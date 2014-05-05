#!/usr/local/bin/perl

#Open the "Non-founder matrix"
open (NFM, "GAW_nonfounder_matrix.txt") || die;


###Load the file into memory

##Reading in non-founder matrix
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
	$NFMlinecount++
}
close (NFM);

###Identify all of the founding chromosomes
my @uniquefounderchrs;
my %seenfounderchrs;

foreach my $value (@NFMfoundingchr) {
	if (! $seenfounderchrs{$value}++ ) {
		push @uniquefounderchrs, $value;
	}
}
print STDOUT "The founding chromosomes are:\n";
foreach $list (@uniquefounderchrs){
print STDOUT "$list\n";
}

###Identify all of the non-founders
my @uniqueindivs;
my %seenindivs;

foreach my $value2 (@NFMindiv) {
	if (! $seenindivs{$value2}++ ) {
		push @uniqueindivs, $value2;
	}
}
print STDOUT "The non-founders are:\n";
foreach $list2 (@uniqueindivs){
print STDOUT "$list2\n";
}

$numuniqueindivs = scalar(@uniqueindivs);
$NFMarraylength = scalar(@NFMindiv);

foreach $currentfoundingchromosome (@uniquefounderchrs){

	###Get just the segments for this foundingchr
	@currentNFMindiv = [];
	@currentNFMfoundingchr = [];
	@currentNFMstart = [];
	@currentNFMend = [];
	@currentNFMPath = [];	
	$currentNFMlinecount = 0;
	for($i = 0; $i < $NFMarraylength; $i++) {
		if (@NFMfoundingchr[$i] eq $currentfoundingchromosome){
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
	

	###Create an array of arrays to store the segment matrix for this chr in
	foreach $NFindiv (@uniqueindivs){
		$currentfoundingchromosome{$NFindiv} = [];
	}
	
	###Identify all possible segment
	my @uniquelocs = ();
	my %seenlocs = ();
	foreach my $loc (@currentNFMstart){
		if(! $seenlocs{$loc}++){
			push @uniquelocs, $loc;
		}
	}
	foreach my $loc2 (@currentNFMend){
		if(! $seenlocs{$loc2}++){
			push @uniquelocs, $loc2;
		}
	}
	
	print STDOUT "Unsorted locs for this NFC:\n";
	foreach $hobo(@uniquelocs){
	print STDOUT "$hobo\n";
	}
	
	my @sorted_locs = ();
	my @sorted_locs = sort {$a <=> $b} @uniquelocs;
	$sortedlocslength = scalar(@sorted_locs);

	print STDOUT "Sorted locs for this NFC:\n";
	foreach $bum(@sorted_locs){
	print STDOUT "$bum\n";
	}

	###Create two arrays with beginning and end of the maximum number of non-overlapping segments
	for($i = 1; $i < $sortedlocslength; $i++) {
		@startarray[$i-1] = @sorted_locs[$i-1];
		@endarray[$i-1] = @sorted_locs[$i];
	}

	###Fill in the array of arrays storing the segment matrix for this chr
	foreach $NFindiv (@uniqueindivs){
		for($i = 0; $i < ($sortedlocslength-1); $i++) {
			my $foundit = 0;
			for($j = 0; $j < $currentNFMarraylength; $j++){
				if(@currentNFMindiv[$j] eq $NFindiv){
					if((@currentNFMstart[$j] < (@startarray[$i]+1))&&(@currentNFMend[$j] > (@endarray[$i]-1))){
						$currentfoundingchromosome{$NFindiv}[$i] = 1;
						$foundit = 1;
					}
				}
			}
			if($foundit == 0){
				$currentfoundingchromosome{$NFindiv}[$i] = 0;
			}
		}
	}
				
	$outfilename = "$currentfoundingchromosome"."_recombination_free_segments.txt";
	open (OUTFILE, ">$outfilename");
	print OUTFILE "Start\tEnd";
	for($i = 0;$i<$numuniqueindivs;$i++){
		print OUTFILE "\t@uniqueindivs[$i]";
	}
	print OUTFILE "\n";
	for($j = 0;$j< ($sortedlocslength-1); $j++){
		print OUTFILE "@startarray[$j]\t@endarray[$j]";
		for($i = 0; $i<$numuniqueindivs;$i++){
			$dude = @uniqueindivs[$i];
			print OUTFILE "\t$currentfoundingchromosome{$dude}[$j]";
		}
		print OUTFILE "\n";
	}

	print STDOUT "The current founding chromosome is $currentfoundingchromosome.\n";
	for($i = 0; $i < ($sortedlocslength-1); $i++) {
		print "@startarray[$i]\t@endarray[$i]\n";
	}

}



