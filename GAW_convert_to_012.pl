#!/usr/local/bin/perl

###Converts ACTG values for genotype data as follows:
###First Homozygous Genotype = 0, Heterozygous Genotype = 1, Second Homozygous Genotype = 2


###ARGV[0] is the genotype data file
open (INFILE, $ARGV[0]) || die;

#open output file converted to 0,1,2's
$filename = $ARGV[0] . "_012.txt";
open (OUTFILE, ">$filename") || die;

#open output file that is a legend of which nucleotide is represented "0" in the 
#"First Homozygous Allele" case
$filename_legend = $ARGV[0] . "_legend.txt";
open (OUTFILE2, ">$filename_legend") || die;

#open output file that contains homozygous SNPs that will be discarded
$filename_homozyg = $ARGV[0] . "_homozyg.txt";
open (OUTFILE3, ">$filename_homozyg") || die;

#read in header line
$line = <INFILE>;
$line =~ s/,/\t/g;
print OUTFILE $line;

chop $line;
@toks = split('\t', $line);
print OUTFILE2 "$toks[0]\t$toks[1]\t$toks[2]\t0\t1\n";

while($line = <INFILE>){
	chop $line;
	@toks = split(',', $line);

	#boolean flag to only write homozygous SNPs to the homozygous file, not the 0,1,2 file
	$homozyg = 1;
	
	#SNP name, Chromosome, and Position to the 012 output array
	$outary[0] = "$toks[0]\t";
	$outary[1] = "$toks[1]\t";
	$outary[2] = "$toks[2]\t";
	#index count to add each call to the array
	$ind = 3;		

	#Choose the first nucleotide in this row to be the "First" for "First Homozygous Allele" as $char0a
	$char0a = substr($toks[3],0,1);

	#For now, assign the same nucleotide to be the "Second" for "Second Homozygous Allele" in case 
	#they are all homozygous with the first nucleotide
	$char0b = substr($toks[3],0,1);
	
	$linelen = @toks;

	for($i = 3; $i < $linelen; $i++){
		$char1 = substr($toks[$i],0,1);
		$char2 = substr($toks[$i],1,1);


		if(($char1 eq $char0a) && ($char2 eq $char0a)){			
			#First Homozygous Genotypes case
			if($char1 eq "-"){
				$outary[$ind] = "-\t";
			}
			else{				
				$outary[$ind] = "0\t";
				#if "--" was the first call for this SNP & there is now a call, adjust for accuracy
				if($char0a eq "-"){
					$char0a = $char1;
					$char0b = $char1;
				}
			}
			$ind++;
		}
		elsif($char1 ne $char2){
			#if "--" was the first call for this SNP & there is now a call, adjust for accuracy
			if($char0a eq "-"){
				$char0a = $char1;
				$char0b = $char2;				
				$outary[$ind] = "1\t";
			}
			else{
				$outary[$ind] = "1\t";
				
				#assign the nucleotide that is the "Second" for "Second Homozygous Allele"
				if($char0a ne $char1){
					$char0b = $char1;
				}
				else{
					 $char0b = $char2;
				}
			}
			$ind++;
			#not homozygous, will output to 0,1,2 file
			$homozyg = 0;
		}
		else{
			#if "--" was the first call for this SNP & there is now a call, adjust for accuracy
			if($char0a eq "-"){
				$char0a = $char1;
				$char0b = $char2;
				$outary[$ind] = "0\t";
			}
			#"--" came after a valid call 
			elsif($char1 eq "-"){
				$outary[$ind] = "-\t";
				$homozyg = 0;
			}
			else{
				#Second Homozygous Allele case
				$outary[$ind] = "2\t";

				#assign the nucleotide that is the "Second" for "Second Homozygous Allele"
				if($char0a ne $char1){
					$char0b = $char1;
				}
				else{
					 $char0b = $char2;
				}
				#not homozygous, will output to 0,1,2 file
				$homozyg = 0;
			}
			$ind++;

		}		
	}

	$outary[$ind] = "\n";
	
	#print out SNP name, Chromosome, and Position to the legend outfile
	print OUTFILE2 "$toks[0]\t$toks[1]\t$toks[2]\t";

	#print out the "First" and "Second" alleles to the legend outfile
	print OUTFILE2 "$char0a\t$char0b\n";

	#if SNP is not homozygous, print the converted data to the output file
	if($homozyg == 0){
		print OUTFILE @outary;
	}
	else{
		print OUTFILE3 "$outary[0]\n";
	}

	#reset homozyg boolean flag for next SNP
	$homozyg = 1;
}
	
	
close INFILE;
close OUTFILE;
close OUTFILE2;
close OUTFILE3;