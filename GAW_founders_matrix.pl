#!/usr/local/bin/perl

### program finds founders, stores the triofile of these founders

#ARGV[0] is pedigree number
$ped = $ARGV[0];

#open list of all trios that have the offspring phased, store in @trios_ary
$filename_cleared = "Cleared_Ped_" . "$ped" . ".txt";
open (INFILE_CLEARED, $filename_cleared) || die;

#open founder outfile
open (OUTFILE, ">GAW_founder_geno.txt") || die;

@trios_ary = <INFILE_CLEARED>;
$filelen = @trios_ary;

###traverse @trios_ary and store triofile for each founder (parent that doesn't have its own parent)

%founder_file_hash = ();
%offspr_hash = ();
@tmp_ary = ();

#store all offspring in %offspr_hash, if is an offspring it is not a founder
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];

	$offspr_hash{$offspr} = 1;
}

#traverse @trios_ary and store file for each founder in %founder_file_hash
#maternal
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];

	if(exists $offspr_hash{$maternal}){
		next;
	}
	else{
		$fname = "$maternal" . "_phased.txt";
		@tmp_ary[0] = $fname;
		@tmp_ary[1] = $maternal;
		@tmp_ary[2] = "m";
		@tmp_ary[3] = "phased";
		if(-e $fname){
			$founder_file_hash{$maternal} = [@tmp_ary];
		}
		else{ 
			$fname = "$maternal" . "_phased_2kids.txt";
			@tmp_ary[0] = $fname;
			@tmp_ary[1] = $maternal;
			@tmp_ary[2] = "m";
			@tmp_ary[3] = "phased";
			if(-e $fname){
				$founder_file_hash{$maternal} = [@tmp_ary];
			}
			else{ 
				$fname = "$trioline" . "_geno.txt";
				@tmp_ary[0] = $fname;
				@tmp_ary[1] = $maternal;
				@tmp_ary[2] = "m";
				@tmp_ary[3] = "m_geno";
				$founder_file_hash{$maternal} = [@tmp_ary];
			}
		}
	}
}

#paternal
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];

	if(exists $offspr_hash{$paternal}){
		next;
	}
	else{
		$fname = "$paternal" . "_phased.txt";
		@tmp_ary[0] = $fname;
		@tmp_ary[1] = $paternal;
		@tmp_ary[2] = "p";
		@tmp_ary[3] = "phased";
		if(-e $fname){
			$founder_file_hash{$paternal} = [@tmp_ary];
		}
		else{
			$fname = "$paternal" . "_phased_2kids.txt";
			@tmp_ary[0] = $fname;
			@tmp_ary[1] = $paternal;
			@tmp_ary[2] = "p";
			@tmp_ary[3] = "phased";
			if(-e $fname){
				$founder_file_hash{$paternal} = [@tmp_ary];
			}
			else{
				$fname = "$trioline" . "_geno.txt";
				@tmp_ary[0] = $fname;
				@tmp_ary[1] = $paternal;
				@tmp_ary[2] = "p";
				@tmp_ary[3] = "p_geno";
				$founder_file_hash{$paternal} = [@tmp_ary];
			}
		}
	}
}

@founder_ary = keys %founder_file_hash;
$num_founders = @founder_ary;

#array for genotype file names for each founder
%fh_ary = ();

#print out column headers
print OUTFILE "SNP\tChr\tPosition\t";

for($i = 0; $i < $num_founders; $i++){	
	@tmp_ary = @{$founder_file_hash{$founder_ary[$i]}};
	@toks = split("_", @tmp_ary[0]);	
	if($tmp_ary[3] eq "m_geno"){
		$id1 = "$toks[0]" . "_" . "$toks[2]" . "_inherited";
		$id2 = "$toks[0]" . "_" . "$toks[2]" . "_not_inherited";
	}
	elsif($tmp_ary[3] eq "p_geno"){ 
		$id1 = "$toks[1]" . "_" . "$toks[2]" . "_inherited";
		$id2 = "$toks[1]" . "_" . "$toks[2]" . "_not_inherited";
	}
	else{
		$id1 = "$toks[0]" . "_1";
		$id2 = "$toks[0]" . "_2";
	}
	print OUTFILE "$id1\t$id2\t";	
	open ($f_infile[$i], $tmp_ary[0]); 
}
print OUTFILE "\n";



#find how many lines are in file
$fname = $f_infile[0];
@lines = <$fname>;
$len = @lines; 
close $f_infile[0];
@tmp_ary = @{$founder_file_hash{$founder_ary[0]}};
open ($f_infile[0], $tmp_ary[0]); 

#read in header line
for($i = 0; $i < $num_founders; $i++){	
	$fname = $f_infile[$i];
	$line = <$fname>;
}

#print out genotypes for founders
for($j = 0; $j < $len; $j++){
	for($i = 0; $i < $num_founders; $i++){	
		$fname = $f_infile[$i];
		$line = <$fname>;
		chop $line;
		@toks = split("\t", $line);

		if($i == 0){
			print OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
		}
	
		@tmp_ary = @{$founder_file_hash{$founder_ary[$i]}};		
		if($tmp_ary[3] eq "m_geno"){
			print OUTFILE "$toks[6]\t$toks[7]\t";
		}
		elsif($tmp_ary[3] eq "p_geno"){
			print OUTFILE "$toks[8]\t$toks[9]\t";
		}
		else{
			print OUTFILE "$toks[3]\t$toks[4]\t";
		}
	}
	print OUTFILE "\n";
}

close OUTFILE;

	


		
































