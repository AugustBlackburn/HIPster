#!/usr/local/bin/perl

#Arguments are the pedigree number and the matrix of diploid genotypes
$ped = $ARGV[0];
$SEQMATRIX = $ARGV[1];

#Open the sequencing data and phased founder chromosomes files
open (SEQdata, $SEQMATRIX) || die;
open (FOUNDERdata, "<Imputed_founder_chrs.txt") || die;

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
@sequencedindivids = @toks2[3 .. $#toks2];

##Read in founder data
print STDOUT "Reading in founder data.\n";

#Get header for founder data and make arrays to hold data
$FOUNDERdataheader = <FOUNDERdata>;
chop $FOUNDERdataheader;
@toks4 = split('\t', $FOUNDERdataheader);
$FOUNDERdataarraylength = scalar (@toks4);
#print STDOUT "Seq data array length is $FOUNDERdataarraylength.\n";
#make array for each column
foreach $colname (@toks4){
	$currentarray{$colname} = [];
	
}
@founderchrsalreadyimputed = @toks4[3 .. $#toks4];

###Identify all of the founders and create an array to hold the founding chromosome data
my @uniquefounder;
my %seenfounder;
@firstfoundingchrs = @toks4[3 .. $#toks4];
foreach my $value (@firstfoundingchrs) {
	@toks6 = split('_',$value);
	$value2 = @toks6[0];
	
	if (! $seenfounder{$value2}++ ) {
		push @uniquefounder, $value2;
	}
}
print STDOUT "The founders are:\n";
foreach $list (@uniquefounder){
print STDOUT "$list\n";
}


foreach $colname (@uniquefounder){
	$founderchr1="$colname"."_1";
	$founderchr2="$colname"."_2";
	$Imputedfounder{$founderchr1} = [];
	$Imputedfounder{$founderchr2} = [];
}

print STDOUT "The founding chromosomes are:\n";
foreach $list (@uniquefounder){
	$founderchr1="$list"."_1";
	$founderchr2="$list"."_2";
	print STDOUT "$founderchr1\n$founderchr2\n";
}

open (FOUNDERinferred, ">Imputed_founder_chrs_2.txt") || die;

print FOUNDERinferred "@toks4[0]\t@toks4[1]\t@toks4[2]";
foreach $colname (@uniquefounder){
	$founderchr1="$colname"."_1";
	$founderchr2="$colname"."_2";
	print FOUNDERinferred "\t$founderchr1\t$founderchr2";
}
print FOUNDERinferred "\n";

while ($SEQdataline = <SEQdata>){
	$FOUNDERdataline = <FOUNDERdata>;
	chop $SEQdataline;
	@toks3 = split('\t', $SEQdataline);
	$SEQdatacolcount = 0;
	foreach $column (@toks3){
		$currentSEQarrayarray = @toks2[$SEQdatacolcount];
		$currentSEQarray{$currentSEQarrayarray}[0]= $column;
		$SEQdatacolcount++;
	}
	$SEQdatalinecount++;
	
	chop $FOUNDERdataline;
	@toks5 = split('\t', $FOUNDERdataline);
	$FOUNDERdatacolcount = 0;
	foreach $column2 (@toks5){
		$currentFOUNDERarrayarray = @toks4[$FOUNDERdatacolcount];
		$currentFOUNDERarray{$currentFOUNDERarrayarray}[0]= $column2;
		$FOUNDERdatacolcount++;
	}
	$FOUNDERdatalinecount++;
	
	foreach $colname (@uniquefounder){
		$founderchr1="$colname"."_1";
		$founderchr2="$colname"."_2";
		$Imputedfounder{$founderchr1}[0] = ();
		$Imputedfounder{$founderchr2}[0] = ();
	}


	foreach $founder (@uniquefounder){
	
	###Is the current founder sequenced?
		$sequenced = 0;
		if (grep {$_ eq $founder} @sequencedindivids){
			$sequenced = 1;
		}
		
		$founderchr1="$founder"."_1";
		$founderchr2="$founder"."_2";
		
		$founderchr1imputed = 0;
		if (grep {$_ eq $founderchr1} @founderchrsalreadyimputed){
			$founderchr1imputed = 1;
		}
		
		$founderchr2imputed = 0;
		if (grep {$_ eq $founderchr2} @founderchrsalreadyimputed){
			$founderchr2imputed = 1;
		}
		
		#print STDOUT "$founder:\t$sequenced\t$founderchr1imputed\t$founderchr2imputed\n";
	
		###for each marker infer allele on paired haplotype if missing, alternatively, if imputed genotypes do not match diploid genotype blank genotypes
		
		
		if($sequenced eq '0'){
			#print STDOUT "Entered founder not sequenced\n";
			if(($founderchr1imputed eq '1') && ($founderchr2imputed eq '1')){
				$Imputedfounder{$founderchr1}[0] = $currentFOUNDERarray{$founderchr1}[0];
				$Imputedfounder{$founderchr2}[0] = $currentFOUNDERarray{$founderchr2}[0];
			}
			elsif(($founderchr1imputed eq '1') && ($founderchr2imputed eq '0')){
				$Imputedfounder{$founderchr1}[0] = $currentFOUNDERarray{$founderchr1}[0];
				$Imputedfounder{$founderchr2}[0] = "u";
			}
			elsif(($founderchr1imputed eq '0') && ($founderchr2imputed eq '1')){
				$Imputedfounder{$founderchr1}[0] = "u";
				$Imputedfounder{$founderchr2}[0] = $currentFOUNDERarray{$founderchr2}[0];
			}
		}
		elsif($sequenced eq '1'){
			if(($founderchr1imputed eq '1') && ($founderchr2imputed eq '1')){
				#print STDOUT "Entered inference 111\n";
				
				$diploid = $currentSEQarray{$founder}[0];
				$founderhaploid1 = $currentFOUNDERarray{$founderchr1}[0];
				$founderhaploid2 = $currentFOUNDERarray{$founderchr2}[0];
				#print STDOUT "$diploid\t$founderhaploid1\t$founderhaploid2\n";
				if($diploid =~ /[0-9]/){
					if(($founderhaploid1 =~ /[0-9]/) && ($founderhaploid2 =~ /[0-9]/)){
						#print STDOUT "All available\n";
						$inferreddiploid = $founderhaploid1 + $founderhaploid2;
						if($diploid == $inferreddiploid){
							$Imputedfounder{$founderchr1}[0] = $founderhaploid1;
							$Imputedfounder{$founderchr2}[0] = $founderhaploid2;
						}
						else{
							$Imputedfounder{$founderchr1}[0] = "u";
							$Imputedfounder{$founderchr2}[0] = "u";
						}
					}
					elsif($founderhaploid1 =~ /[0-9]/){
						#print STDOUT "founder2 not available\n";
						$founderhaploid2 = $diploid - $founderhaploid1;
						$Imputedfounder{$founderchr1}[0] = $founderhaploid1;
						$Imputedfounder{$founderchr2}[0] = $founderhaploid2;
					}
					elsif($founderhaploid2 =~ /[0-9]/){
						#print STDOUT "founder1 not available\n";
						$founderhaploid1 = $diploid - $founderhaploid2;
						$Imputedfounder{$founderchr1}[0] = $founderhaploid1;
						$Imputedfounder{$founderchr2}[0] = $founderhaploid2;
					}
					else{
						#print STDOUT "Neither founder available\n";
						$Imputedfounder{$founderchr1}[0] = $currentFOUNDERarray{$founderchr1}[0];
						$Imputedfounder{$founderchr2}[0] = $currentFOUNDERarray{$founderchr2}[0];
					}
				}
				else{
					#print STDOUT "Diploid not available\n";
					$Imputedfounder{$founderchr1}[0] = $founderhaploid1;
					$Imputedfounder{$founderchr2}[0] = $founderhaploid2;
				}
	
			}
			elsif(($founderchr1imputed eq '1') && ($founderchr2imputed eq '0')){
				$diploid = $currentSEQarray{$founder}[0];
				$founderhaploid1 = $currentFOUNDERarray{$founderchr1}[0];
							
				if($diploid =~ /[0-9]/){
				#print STDOUT "Entered inference 110\n";
					if($founderhaploid1 =~ /[0-9]/){
						$founderhaploid2 = $diploid - $founderhaploid1;	
						$Imputedfounder{$founderchr1}[0] = $founderhaploid1;
						$Imputedfounder{$founderchr2}[0] = $founderhaploid2;
						
					}
					else{
						$Imputedfounder{$founderchr1}[0] = "u";
						$Imputedfounder{$founderchr2}[0] = "u";
					}
				}
				else{
					$Imputedfounder{$founderchr1}[0] = $currentFOUNDERarray{$founderchr1}[0];
					$Imputedfounder{$founderchr2}[0] = "u";
				}			
			}
			elsif(($founderchr1imputed eq '0') && ($founderchr2imputed eq '1')){
				#print STDOUT "Entered inference 110\n";
				$diploid = $currentSEQarray{$founder}[0];
				$founderhaploid2 = $currentFOUNDERarray{$founderchr2}[0];
				$founderhaploid1 = $diploid - $founderhaploid2;
				$Imputedfounder{$founderchr1}[0] = $founderhaploid1;
				$Imputedfounder{$founderchr2}[0] = $founderhaploid2;	
				
			}
		}	
	}
	print FOUNDERinferred "$currentFOUNDERarray{@toks4[0]}[0]\t$currentFOUNDERarray{@toks4[1]}[0]\t$currentFOUNDERarray{@toks4[2]}[0]";
	foreach $colname (@uniquefounder){
	$founderchr1="$colname"."_1";
	$founderchr2="$colname"."_2";
	print FOUNDERinferred "\t$Imputedfounder{$founderchr1}[0]";
	print FOUNDERinferred "\t$Imputedfounder{$founderchr2}[0]";
	}
	print FOUNDERinferred "\n";
}
close (FOUNDERinferred);
close (SEQdata);
close (FOUNDERdata);
	