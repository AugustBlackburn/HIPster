#!/usr/local/bin/perl

use Switch;


#$ARGV[0] is the converted 012 file
open (INFILE_GENO, $ARGV[0]) || die;

#$ARGV[1] is pedigree file
open (INFILE_PED, $ARGV[1]) || die;

#$ARGV[2] is the pedigree number
$ped = $ARGV[2];

#file of SNPs with genotype errors
$filename_genoerr = "Genotype_errors" . "$ped" . ".txt";
open (OUTFILE_GENOERR, ">$filename_genoerr") || die;

#cleared trios file
$filename_cleared = "Cleared_Ped_" . "$ped" . ".txt";
open (OUTFILE_CLEARED, ">$filename_cleared") || die;

#Open geno error count file and create header
open (OUTFILE_GENOERR_COUNT, ">Genotype_error_count.txt") || die;
print OUTFILE_GENOERR_COUNT "Mother\tFather\tID\tGeno_error_count\tHeterozygosity_ambiguity_count\tMissing_count\n";

#make a hash of indices of the individual id's from the genotype header line
%indices = ();
$line = <INFILE_GENO>;
chop $line;
@toks = split('\t', $line);
$ct = 3;
$linelen = @toks;
close INFILE_GENO;

for ($i = 3; $i < $linelen; $i++){
	$indices{$toks[$i]} = $i;	
}

$error_message = "1";
$missing_message = "2";
#traverse Pedigree file for useable trios, read in header line first
$line = <INFILE_PED>;

while($line = <INFILE_PED>){
	chop $line;
	@toks = split(',', $line);
	#check if is a trio	
	if(($toks[1] && $toks[2])){
		$offspr = $toks[0];
		$paternal = $toks[1];
		$maternal = $toks[2];

		$trio_id = "$maternal"."_"."$paternal"."_"."$offspr";
				
		$m_ind = $indices{$maternal};
		$p_ind = $indices{$paternal};
		$o_ind = $indices{$offspr};
		
		#if there is genotype data for the trio, output the trio file, otherwise skip to next line in ped file
		if($m_ind && $p_ind && $o_ind){		
			print OUTFILE_CLEARED "$trio_id\n";
			open (INFILE_GENO, $ARGV[0]);
			#read in header line of genotype file;
			$genoline = <INFILE_GENO>;

			$filename_Trio = "$trio_id" . "_geno.txt";
			open (OUTFILE_TRIO, ">$filename_Trio");
			print OUTFILE_TRIO "SNP Name\tChr\tPosition\t" . "$maternal\t$paternal\t$offspr\t" . "$offspr"."_"."$maternal"."_inherited\t"
					. "$offspr"."_"."$maternal"."_not_inherited\t" . "$offspr"."_"."$paternal"."_inherited\t" . "$offspr"."_"."$paternal"."not_inherited\n";
			$genoerrorcount = 0;
			$heterozygositycount = 0;
			$missingerror = 0;

			while($genoline = <INFILE_GENO>){
				chop $genoline;
				@toks = split('\t', $genoline);

				$trio_str = $toks[$m_ind] . $toks[$p_ind] . $toks[$o_ind];

				if($trio_str =~ m/-/){
					$trio_str = "-";
				}				
				
				switch($trio_str){
					case "000"	{$out_str = "0\t0\t0\t0\t0\t0\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "021"	{$out_str = "0\t2\t1\t0\t0\t1\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "201"	{$out_str = "2\t0\t1\t1\t1\t0\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "222"	{$out_str = "2\t2\t2\t1\t1\t1\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "010"	{$out_str = "0\t1\t0\t0\t0\t0\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "011"	{$out_str = "0\t1\t1\t0\t0\t1\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "100"	{$out_str = "1\t0\t0\t0\t1\t0\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "101"	{$out_str = "1\t0\t1\t1\t0\t0\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "211"	{$out_str = "2\t1\t1\t1\t1\t0\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "212"	{$out_str = "2\t1\t2\t1\t1\t1\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "121"	{$out_str = "1\t2\t1\t0\t1\t1\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "122"	{$out_str = "1\t2\t2\t1\t0\t1\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "110"	{$out_str = "1\t1\t0\t0\t1\t0\t1\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "112"	{$out_str = "1\t1\t2\t1\t0\t1\t0\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;}

					case "111"	{$out_str = "1\t1\t1\tu\tu\tu\tu\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;
							print OUTFILE_GENOERR "$toks[0]\t$toks[1]\t$toks[2]\t$trio_str\t0\t$maternal\t$paternal\t$offspr\n";
							$heterozygositycount++ ;}

					case "-"	{$out_str = "m\tm\tm\tm\tm\tm\tm\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str;
							print OUTFILE_GENOERR "$toks[0]\t$toks[1]\t$toks[2]\t$trio_str\t$missing_message\t$maternal\t$paternal\t$offspr\n"; 
							$missingerror++ ;}
			
					else		{$out_str = "e\te\te\te\te\te\te\t\n";
							print OUTFILE_TRIO "$toks[0]\t$toks[1]\t$toks[2]\t";
							print OUTFILE_TRIO $out_str; 
							print OUTFILE_GENOERR "$toks[0]\t$toks[1]\t$toks[2]\t$trio_str\t$error_message\t$maternal\t$paternal\t$offspr\n";
							$genoerrorcount++ ;}						

				}
				
			}
			print OUTFILE_GENOERR_COUNT "$maternal\t$paternal\t$offspr\t$genoerrorcount\t$heterozygositycount\t$missingerror\n";
			close INFILE_GENO;
			close OUTFILE_TRIO;
		}			
	}

}

close INFILE_PED;
close INFILE_GENO;
close OUTFILE_CLEARED;
close OUTFILE_GENOERR;
close OUTFILE_GENOERR_COUNT;

#%bad_snps = ();

#open (INFILE_GENOERR, "$filename_genoerr") || die;

#while($badline = <INFILE_GENOERR>){
#	@toks = split('\t', $badline);
#	$len = @toks;
#	#is a non "111" error
#	if($len == 5){
#		$bad_snps{$toks[0]} = 1;
#	}
#}

#@files = <*_geno.txt>;
#foreach $file (@files){
#	$tempfile = "$file" . "_temp.txt";

#	open (INFILE_GENO, $file) || die;
#	open (OUTFILE_GENO, ">$tempfile") || die;

#	#take care of header line
#	$line = <INFILE_GENO>;	
#	print OUTFILE_GENO $line;
#	while($line = <INFILE_GENO>){
#		@toks = split('\t', $line);
#		if(exists $bad_snps{$toks[0]}){
#			next;
#		}
#		else{
#			print OUTFILE_GENO $line;
#		}
#	}
#	close INFILE_GENO;
#	close OUTFILE_GENO;
#	rename $tempfile, $file;
#}











