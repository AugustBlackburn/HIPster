#!/usr/local/bin/perl

#Arguments are the pedigree number and the matrix of diploid genotypes
$ped = $ARGV[0];
$SEQMATRIX = $ARGV[1];

#Open the "Non-founder matrix", "clearedped", and the sequencing data files
open (NFM, "GAW_nonfounder_matrix.txt") || die;
$filename_cleared = "Cleared_Ped_"."$ped".".txt";
open (CLEAREDPED, $filename_cleared) || die;
open (SEQdata, $SEQMATRIX) || die;

###Load the files into memory

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
	$NFMlinecount++
}
close (NFM);


##Read in seq data
print STDOUT "Reading in SEQ data.\n";

#Get header for seq data and make arrays to hold data
$SEQdataheader = <SEQdata>;
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
while($SEQdataline = <SEQdata>){
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
print "There are $SEQdatalinecount markers in $SEQMATRIX.\n";
$numberofsequencedindividuals = $SEQdataarraylength - 3;
print "There are $numberofsequencedindividuals available individuals with unphased sequence data available.\n";
close (SEQdata);


##Read in the cleared-ped file and split into mother, father, offspring
print STDOUT "Reading in cleared ped file.\n";
$triolinecount = 0;
while($trioline = <CLEAREDPED>){
	chop $trioline;
	@toks4 = split('_',$trioline);
	@mother[$triolinecount] = @toks4[0];
	@father[$triolinecount] = @toks4[1];
	@offspring[$triolinecount] = @toks4[2];
	$triolinecount++;
}
print "There are $triolinecount individuals with phased sequence data available.\n";
close (CLEAREDPED);

####This is code for checking to make sure the data was read in correctly.
$NFMarraylength = scalar(@NFMindiv);
#print STDOUT "NFM array length is $NFMarraylength.\n";
#for($i = 0; $i < $NFMarraylength; $i++){
#	print STDOUT "@NFMindiv[$i],@NFMfoundingchr[$i],@NFMstart[$i],@NFMend[$i],@NFMPath[$i]\n";
#}




#for( $i = 0; $i < 50; $i++) {
#	foreach $colname (@toks2){
#		print "$currentarray{$colname}[$i]\t";
#	}
#}

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

####Create list of all sequenced individuals
@sequencedindivids = @toks2[3 .. $#toks2];
##This is a check to make sure individs are correct
#foreach $dude (@sequencedindivids){
#	print "$dude\n";
#}

####For each founding chromosome, use available information to impute genotypes

###Create an array of arrays to store the phased genotypes for each founding chr
foreach $founderchrforarray (@uniquefounderchrs){
	$arrayoffounderchrs{$founderchrforarray} = [];
	
}
###Set global variable for blanked and uninformative genotypes
$blankgenotype = "-";
$uninformative = "u";

foreach $currentfoundingchromosome (@uniquefounderchrs){
	####Is current founder sequenced
	$sequenced = 0;
	@toks5 = split('_',$currentfoundingchromosome);
	$currfounder = @toks5[0];
	if (grep {$_ eq $currfounder} @sequencedindivids){
		$sequenced = 1;
	}
	
	print STDOUT "The founder ($currfounder) of the current founding chromosome ($currentfoundingchromosome) sequence status is $sequenced.\n";
	print STDOUT "Imputing genotypes for $currentfoundingchromosome. Available sequence data for $currfounder will be used.\n";
	####Create shortened non-founder matrix for currentfoundingchr
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
	#print STDOUT "The current NFMarray length is $currentNFMarraylength.\n";
	###print statement for checking that the shortened NFM is working
	#print "$currentfoundingchromosome is the current founding chromosome.\n\n";
	#for($i = 0; $i < $NFMarraylength; $i++) {
	#	print "@currentNFMindiv[$i]\t@currentNFMfoundingchr[$i]\t@currentNFMstart[$i]\t@currentNFMend[$i]\t@currentNFMPath[$i]\n";
	#}
	#print "\n";
	
	#Identify non-founders who inherited segment of current founder chromosome.
	my @uniquenonfounderchrs = [];
	my %seennonfounderchrs = [];

		foreach $NFvalue (@currentNFMindiv) {
			#print STDOUT "$NFvalue\n";
			if (! $seennonfounderchrs{$NFvalue}++ ) {
				push @uniquenonfounderchrs, $NFvalue;
			}
		}


	#print STDOUT "The non-founders carrying a segment of $currentfoundingchromosome are:\n";
	#foreach $NFlist (@uniquenonfounderchrs){
	#print STDOUT "$NFlist\n";
	#}
	
	#Identify which of these individuals is sequenced
	my @sequencedNFforthisNFchr = [];
	$seqindivtok = 0;
	foreach $NFseg (@uniquenonfounderchrs){
		foreach $seqedindiv (@sequencedindivids) {
			if ($NFseg eq $seqedindiv) {
				@sequencedNFforthisNFchr[$seqindivtok]= $seqedindiv;
				$seqindivtok++;
			}
		}
	}
	#print STDOUT "@sequencedNFforthisNFchr\n";
	
	#Identfy non-founders who are phased
	my @phasedNFforthisNFchr = [];
	$phaseindivtok = 0;
	foreach $NFphase (@sequencedNFforthisNFchr){
		foreach $phasedindiv (@offspring) {
			#print STDOUT "$NFphase might be $phasedindiv\n";
			if ($NFphase eq $phasedindiv) {
				#print STDOUT "$NFphase is $phasedindiv\n";
				@phasedNFforthisNFchr[$phaseindivtok] = $NFphase;
				$phaseindivtok++;
			}
		}
	}
	#print STDOUT "@phasedNFforthisNFchr\n";

	#make array of arrays to hold phasedNFdata
	foreach $currphasekid (@phasedNFforthisNFchr){
	$phasedkidarray{$currphasekid} = [];
	}
	#count number of phased kids available for the 	NF chr
	$phasedkidforthisNFchrcount = scalar(@phasedNFforthisNFchr);
	#print STDOUT "There are $phasedkidforthisNFchrcount phased individuals for this NF chr.\n";

	###Open files for phased non-founders
	foreach $phasedkid (@phasedNFforthisNFchr) {
		#print STDOUT "the current kid is $phasedkid\n";
		for ($i=0; $i<$triolinecount; $i++){
			if(@offspring[$i] eq $phasedkid){
				my $phasedkidfilename = ("@mother[$i]"."_"."@father[$i]"."_"."@offspring[$i]"."_geno.txt");
				#print STDOUT "the file name is $phasedkidfilename\n";
				open (PHASEDKIDFILE, $phasedkidfilename) || die;
				$phasedkidheader = <PHASEDKIDFILE>;
				@toks8 = split('\t',$phasedkidheader);
				#print STDOUT "the current kids header is $phasedkidheader\n\n";
				@currentNFindiv = [];
				@currentNFfoundingchr = [];
				@currentNFstart = [];
				@currentNFend = [];
				@currentNFPath = [];	
				$currentNFlinecount = 0;
				for($i = 0; $i < $currentNFMarraylength; $i++) {
					#print STDOUT "@currentNFMindiv[$i] might be $phasedkid.\n";
					if (@currentNFMindiv[$i] eq $phasedkid){
						#print STDOUT "@currentNFMindiv[$i] IS $phasedkid.\n";
						#print STDOUT "@NFMfoundingchr[$i] = $currentfoundingchromosome\n";
						@currentNFindiv[$currentNFlinecount] = @currentNFMindiv[$i];
						@currentNFfoundingchr[$currentNFlinecount] = @currentNFMfoundingchr[$i];
						@currentNFstart[$currentNFlinecount] = @currentNFMstart[$i];
						@currentNFend[$currentNFlinecount] = @currentNFMend[$i];
						@currentNFPath[$currentNFlinecount] = @currentNFMPath[$i];
						$currentNFlinecount++;
					}
				}
				##Read in the correct inherited chromosome for the kid and store in array
				$currentNFPathtoexplode = @currentNFPath[0];
				@toks6 = split('_',$currentNFPathtoexplode);
				$checkforinherited = @currentNFMfoundingchr[0];
				@toks7 = split('_',$checkforinherited);
				$inheritedword = "inherited";
				if(@toks7[2] eq $inheritedword){
					$colnameofinterest = ("@toks6[1]"."_"."@toks6[0]"."_inherited");
					#print STDOUT "$colnameofinterest\n";
					
				}
				else{
					$colnameofinterest = ("@toks6[0]"."_"."@toks6[1]"."_inherited");
					#print STDOUT "$colnameofinterest\n";
				}
				$k = [];
				for($j = 0; $j <10; $j++){
					$check = @toks8[$j];
					if ($check eq $colnameofinterest){
						$k = $j;
					}
				}
				while ($line = <PHASEDKIDFILE>){
					chop $line;
					@toks9 = split("\t",$line);
					push @{$phasedkidarray{$phasedkid}},@toks9[$k];
				}
				for($i = 0; $i <50; $i++){
				#print STDOUT "$phasedkidarray{$phasedkid}[$i]\n";
				}

				close (PHASEDKIDFILE);
			}
		}
		
	}
	#####For each marker follow decision tree
	#####Decision tree for most likely genotype for this founder chromosome...
		#### Use boolean flag boolean = 0
		#### Use boolean flag boolean2 = 0
		#### Is the founder sequenced?
			###yes
				####is the genotype unambiguous ie. 0 or 2? 
					###yes
						####use genotype
					###no?
						#### boolean = 1
			###no
				#### boolean = 1
		####does boolean = 1
			###yes 
				####is there a phased result?
					###yes
						###use phased result
					###no
						###boolean2 = 1
			###no
				####next
		####does boolean2 = 1
			###yes
				####is there an unambigous genotype from NF sequence ie. 0 or 2 and no conflict?
					###yes
						###use unambiguous result
					###no
						###output "u" for uninformative
			###no
				####next
	$currentmarkercount = 0;
	$locarrayname = @toks2[2];
	foreach $currmarkerloc (@{$currentarray{$locarrayname}}){
		$boolean1 = 0;
		$boolean2 = 0;		
		if ($sequenced == 1){
			$sequencedmarkergeno = $currentarray{$currfounder}[$currentmarkercount];
			#print STDOUT "$currmarkerloc\t$sequencedmarkergeno\n";
			if ($sequencedmarkergeno eq $blankgenotype){
				#print STDOUT "$currmarkerloc is a blank genotype for $currfounder.\n";
				$boolean1 = 1;
			}
			elsif($sequencedmarkergeno eq 1){
				$boolean1 = 1;
			}
			elsif($sequencedmarkergeno eq 0){
				push @{$arrayoffounderchrs{$currentfoundingchromosome}},0;
			}
			elsif($sequencedmarkergeno eq 2){
				push @{$arrayoffounderchrs{$currentfoundingchromosome}},1;
			}
			else{
				$boolean1 = 1;
			}

		}
		else {
			$boolean1 = 1;
		}
		if ($boolean1 == 1){
			#print STDOUT "Entered boolean1 at $currmarkerloc.\n";
			@arrayofobservedphasedgenotypes = [];
			$counter1 = 0;
			foreach $phasedkidholder(@phasedNFforthisNFchr){
				for($i = 0; $i < $currentNFMarraylength; $i++) {
					if(@currentNFMindiv[$i] eq $phasedkidholder){
						if(@currentNFMstart[$i] <= $currmarkerloc){
							if(@currentNFMend[$i] >= $currmarkerloc){
								push @arrayofobservedphasedgenotypes, $phasedkidarray{$phasedkidholder}[$currentmarkercount];
								$counter1++;
							}
						}
					}
				}
			}
			if($counter1 == 0){
				$boolean2 = 1;
			}
			elsif($counter1 > 0){
				#print STDOUT "Determining position $currmarkerloc from phased data.\n";
				$phasezerocount = 0;
				$phaseonecount = 0;
				foreach $phasegenotypetoken (@arrayofobservedphasedgenotypes){
					#print STDOUT "$phasegenotypetoken\n";
					if($phasegenotypetoken eq 0){
						$phasezerocount++;
					}
					elsif($phasegenotypetoken eq 1){
						$phaseonecount++;
					}
				}
				#print STDOUT "Phase zero count:$phasezerocount.\n";
				#print STDOUT "Phase one  count:$phaseonecount.\n";
				if($phasezerocount == 0){
					if($phaseonecount == 0){
						$boolean2 = 1;
					}
					elsif($phaseonecount >= 1){
						push @{$arrayoffounderchrs{$currentfoundingchromosome}},1;
					}
				}
				elsif($phasezerocount > 0){
					if($phaseonecount == 0){
						push @{$arrayoffounderchrs{$currentfoundingchromosome}},0;	
					}
					elsif($phaseonecount >= 1){
						$boolean2 = 1;
					}
				}
					
			}	

		}
		if ($boolean2 == 1){
			#print STDOUT "Entered boolean2 at $currmarkerloc.\n";
			@arrayofobservedsequencedgenotypes = [];
			$counter2 = 0;
			foreach $sequencedkidholder (@sequencedNFforthisNFchr){	
				for($i = 0; $i < $currentNFMarraylength; $i++) {
					if(@currentNFMindiv[$i] eq $sequencedkidholder){	
						if(@currentNFMstart[$i] <= $currmarkerloc){
							if(@currentNFMend[$i] >= $currmarkerloc){
								push @arrayofobservedsequencedgenotypes, $currentarray{$sequencedkidholder}[$currentmarkercount];
								$counter2++;
							}
						}
					}
				}
			}
			if($counter2 ==0){
				push @{$arrayoffounderchrs{$currentfoundingchromosome}},$uninformative;
			}
			elsif($counter2 > 0){
				#print STDOUT "Determining position $currmarkerloc from sequenced non-founder data.\n";
				$sequencedzerocount = 0;
				$sequencedtwocount = 0;
				foreach $sequencedgenotypetoken (@arrayofobservedsequencedgenotypes){
					#print STDOUT "$sequencedgenotypetoken";
					if($sequencedgenotypetoken eq 0){
						$sequencedzerocount++;
					}
					elsif($sequencedgenotypetoken eq 2){
						$sequencedtwocount++;
					}
				}
				#print STDOUT "\nSeq zero count:$sequencedzerocount.\n";
				#print STDOUT "Seq two count:$sequencedtwocount.\n";
				if($sequencedzerocount == 0){
					if($sequencedtwocount == 0){
						push @{$arrayoffounderchrs{$currentfoundingchromosome}},$uninformative;
					}
					elsif($sequencedtwocount >=1){
						push @{$arrayoffounderchrs{$currentfoundingchromosome}},1;
					}
				}
				elsif($sequencedzerocount > 0){
					if($sequencedtwocount == 0){
						push @{$arrayoffounderchrs{$currentfoundingchromosome}},0;
					}
					elsif($sequencedtwocount >=1){
						push @{$arrayoffounderchrs{$currentfoundingchromosome}},$uninformative;
					}
				}
			}

		}
	$currentmarkercount++;
	}
}
for ($j=3;$j<$SEQdataarraylength;$j++){
	$colheader3 = @toks2[$j];
	print STDOUT "$currentarray{$colheader3}\n";
	delete $currentarray{$colheader3};
}

undef @NFMindiv;
undef @NFMfoundingchr;
undef @NFMstart;
undef @NFMend;
undef @NFMPath;
undef @currentNFMindiv;
undef @currentNFMfoundingchr;
undef @currentNFMstart;
undef @currentNFMend;
undef @currentNFMPath;

print STDOUT "Finished imputing into founding chromosomes. Writing results to Imputed_founder_chrs.txt.\n";
#output phased founder chrs
open IMPUTEDFOUNDERCHR, ">Imputed_founder_chrs.txt" or die $!;
foreach $colheader (@toks2[0..2]){
	print IMPUTEDFOUNDERCHR "$colheader\t";
}
$colheader1 = @toks2[0];
$colheader2 = @toks2[1];
$colheader3 = @toks2[2];
foreach $outputfounderchr (@uniquefounderchrs){
	print IMPUTEDFOUNDERCHR "$outputfounderchr\t";
}
print IMPUTEDFOUNDERCHR "\n";
#print STDOUT "There are $SEQdatalinecount markers to output.";
for($i=0;$i < $SEQdatalinecount; $i++){
	print IMPUTEDFOUNDERCHR "$currentarray{$colheader1}[$i]\t$currentarray{$colheader2}[$i]\t$currentarray{$colheader3}[$i]\t";
	foreach $outputfounderchr (@uniquefounderchrs){
		print IMPUTEDFOUNDERCHR "$arrayoffounderchrs{$outputfounderchr}[$i]\t";
		delete $arrayoffounderchrs{$outputfounderchr}[$i];
	}
	print IMPUTEDFOUNDERCHR "\n";
}
close (IMPUTEDFOUNDERCHR);