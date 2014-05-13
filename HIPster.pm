sub Step1 {
#@_[0] is the converted 012 file, @_[1] is pedigree file
open (INFILE_GENO, @_[0]) || die;
open (INFILE_PED, @_[1]) || die;

open (OUTFILE_GENOERR, ">Genotype_errors.csv") || die;
open (OUTFILE_CLEARED, ">Cleared_trios.csv") || die;
open (OUTFILE_GENOERR_COUNT, ">Genotype_error_count.csv") || die;
print OUTFILE_GENOERR_COUNT "Mother,Father,ID,Geno_error_count,Heterozygosity_ambiguity_count,Missing_count\n";

#make a hash of indices of the individual id's from the genotype header line
%indices = ();
$line = <INFILE_GENO>;
chop $line;
@toks1 = split(',', $line);
$length = scalar(@toks1);
close INFILE_GENO;

for ($i = 3; $i < $length; $i++){
	$indices{$toks1[$i]} = $i;	
}

$error_message = "1";
$missing_message = "2";
#traverse Pedigree file for useable trios, read in header line first
$headerline = <INFILE_PED>;

while($line = <INFILE_PED>){
	chop $line;
	@toks2 = split(',', $line);
	#check if is a trio	
	if(($toks2[1] && $toks2[2])){
		$offspr = $toks2[0];
		$paternal = $toks2[1];
		$maternal = $toks2[2];

		$trio_id = "$maternal"."_"."$paternal"."_"."$offspr";
				
		$m_ind = $indices{$maternal};
		$p_ind = $indices{$paternal};
		$o_ind = $indices{$offspr};
		
		#if there is genotype data for the trio, output the trio file, otherwise skip to next line in ped file
		if($m_ind && $p_ind && $o_ind){
			print OUTFILE_CLEARED "$trio_id\n";
			open (INFILE_GENO, @_[0]);
			#read in header line of genotype file;
			$genoline = <INFILE_GENO>;

			$filename_Trio = "$trio_id" . "_geno.csv";
			open (OUTFILE_TRIO, ">$filename_Trio");
			print OUTFILE_TRIO "SNP Name,Chr,Position," . "$maternal,$paternal,$offspr," . "$offspr"."_"."$maternal"."_inherited,"
					. "$offspr"."_"."$maternal"."_not_inherited," . "$offspr"."_"."$paternal"."_inherited," . "$offspr"."_"."$paternal"."not_inherited\n";
			$genoerrorcount = 0;
			$heterozygositycount = 0;
			$missingerror = 0;

			while($genoline = <INFILE_GENO>){
				chop $genoline;
				@toks3 = split(',', $genoline);

				$trio_str = $toks3[$m_ind] . $toks3[$p_ind] . $toks3[$o_ind];

				if($trio_str =~ m/-/){
					$trio_str = "-";
				}				
				
				if($trio_str eq "000"){
					$out_str = "0,0,0,0,0,0,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "021"){
					$out_str = "0,2,1,0,0,1,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "201"){
					$out_str = "2,0,1,1,1,0,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "222"){
					$out_str = "2,2,2,1,1,1,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "010"){
					$out_str = "0,1,0,0,0,0,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "011"){
					$out_str = "0,1,1,0,0,1,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "100"){
					$out_str = "1,0,0,0,1,0,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "101"){
					$out_str = "1,0,1,1,0,0,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "211"){
					$out_str = "2,1,1,1,1,0,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "212"){
					$out_str = "2,1,2,1,1,1,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "121"){
					$out_str = "1,2,1,0,1,1,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "122"){
					$out_str = "1,2,2,1,0,1,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "110"){
					$out_str = "1,1,0,0,1,0,1\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "112"){
					$out_str = "1,1,2,1,0,1,0\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
				}elsif($trio_str eq "111"){
					$out_str = "1,1,1,u,u,u,u\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
					print OUTFILE_GENOERR "$toks3[0],$toks3[1],$toks3[2],$trio_str,0,$maternal,$paternal,$offspr\n";
					$heterozygositycount++ ;
				}elsif($trio_str eq "-"){
					$out_str = "m,m,m,m,m,m,m\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str;
					print OUTFILE_GENOERR "$toks3[0],$toks3[1],$toks3[2],$trio_str,$missing_message,$maternal,$paternal,$offspr\n"; 
					$missingerror++ ;
				}else{
					$out_str = "e,e,e,e,e,e,e\n";
					print OUTFILE_TRIO "$toks3[0],$toks3[1],$toks3[2],";
					print OUTFILE_TRIO $out_str; 
					print OUTFILE_GENOERR "$toks3[0],$toks3[1],$toks3[2],$trio_str,$error_message,$maternal,$paternal,$offspr\n";
					$genoerrorcount++ ;
				}										
			}
			print OUTFILE_GENOERR_COUNT "$maternal,$paternal,$offspr,$genoerrorcount,$heterozygositycount,$missingerror\n";
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
}

sub Step2 {
#@_[0] is the number of markers required for a switch to be considered "real"

print STDOUT "The number of markers required for a switch is @_[0].\n";
if(@_[0] == 0){
	print STDOUT "Please enter a value greater than 0 for the number of markers necessary for switches to be considered real\n";
	die;
}
if(@_[0] > 10){
	print STDOUT "You have indicated that more than 10 markers are necessary to support a recombination event. Consider using less markers if you obtain odd results.\n";
}

@files = glob("*3gen.csv"); 
unlink @files;

#open list of all trios that have the offspring phased, store in @trios_ary
open (INFILE_CLEARED, "Cleared_trios.csv") || die;

###Read in cleared trios file and store hashes with mat and pat as keys, trio file(s) as value(s) if the mat or pat is also an offspring
%offspr_hash = ();
%offspr_mat_hash = ();
%offspr_pat_hash = ();

#first collect all offspr into $offspr_hash
while($trioline=<INFILE_CLEARED>){	
	chomp $trioline;
	@toks1 = split('_', $trioline);
	$maternal = @toks1[0];
	$paternal = @toks1[1];
	$offspr = @toks1[2];
	push @mats, $maternal;
	push @pats, $paternal;
	push @kids, $offspr;
	@{$offspr_hash{$offspr}}[0] = $trioline;
	push @{$offspr_pat_hash{$paternal}}, $offspr;
	push @{$offspr_mat_hash{$maternal}}, $offspr;
}

#write 3_gen_crossover file for paternal
foreach $offspr (keys %offspr_hash){
	if(exists $offspr_pat_hash{$offspr}){
		$numkids = scalar(@{$offspr_pat_hash{$offspr}});
		#for each offspring of a father who is also a kid in another trio, open the fathers trio file and the offsprings trio file
		@toks=split('_',@{$offspr_hash{$offspr}}[0]);
		$gm_pat = @toks[0];
		$gp_pat = @toks[1];
		foreach $offspring (@{$offspr_pat_hash{$offspr}}){
			$fname_out = "$offspring" . "_crossovers_3gen.csv";
			open (OUTFILE, ">$fname_out") || die;
			$offspring_trio_filename = "@{$offspr_hash{$offspring}}[0]"."_geno.csv";
			$paternal_trio_filename = "@{$offspr_hash{$offspr}}[0]"."_geno.csv";
			open (INFILE_P, "$paternal_trio_filename") || die;
			open (INFILE_O, "$offspring_trio_filename") || die;
			#read in the headers for each trio file
			$headerline1=<INFILE_P>;
			$headerline2=<INFILE_O>;
			while($linep=<INFILE_P>){
				$lineo=<INFILE_O>;
				chomp $linep;
				chomp $lineo;
				@toksp = split(',',$linep);
				@tokso = split(',',$lineo);
				if((join(',',@toksp[3..9]) =~ /[emu]/) || (join(',',@tokso[3..9]) =~ /[emu]/)){
					next;
				}
				if(@toksp[5] == 1){
					if(@toksp[8] == @tokso[8]){
						print OUTFILE "$gp_pat"."_"."$offspr,";
						print OUTFILE "@tokso[2],";
						$currg = $gp_pat;
						$prev_loc = @tokso[2];
						last;
					}elsif(@toksp[6] == @tokso[8]){
						print OUTFILE "$gm_pat" . "_" ."$offspr,";
						print OUTFILE "@tokso[2],";
						$currg = $gm_pat;
						$prev_loc = @tokso[2];
						last;
					}
				}
			}#found first informative SNP and exited while loop
			$switch_now = 0;
			while($linep=<INFILE_P>){
				$lineo=<INFILE_O>;
				chomp $linep;
				chomp $lineo;
				@toksp = split(',',$linep);
				@tokso = split(',',$lineo);
				if((join(',',@toksp[3..9]) =~ /[emu]/) || (join(',',@tokso[3..9]) =~ /[emu]/)){
					next;
				}
				if(@toksp[5] == 1){
					if(@toksp[6] == @tokso[8]){
						if($currg eq $gm_pat){
							$prev_loc = @tokso[2];
							$switch_now = 0;
							next;
						}else{
							if($switch_now == 0){
								$end_loc = $prev_loc;
								$start_loc = @tokso[2];
							}
							$switch_now++;
							$prev_loc = @tokso[2];
							if($switch_now==@_[0]){
								print OUTFILE "$end_loc\n$gm_pat"."_"."$offspr,$start_loc,";
								$currg = $gm_pat;
								$switch_now=0;
							}
						}
					}elsif(@toksp[8] == @tokso[8]){
						if($currg eq $gp_pat){
							$prev_loc=@tokso[2];
							$switch_now = 0;
							next;
						}else{
							if($switch_now == 0){
								$end_loc=$prev_loc;
								$start_loc=@tokso[2];
							}
							$switch_now++;
							$prev_loc = @tokso[2];
							if($switch_now == @_[0]){
								print OUTFILE "$end_loc\n$gp_pat"."_"."$offspr,$start_loc,";
								$currg = $gp_pat;
								$switch_now = 0;
							}
						}
					}
				}
			}
			print OUTFILE "$prev_loc\n";
			close INFILE_O;
			close INFILE_P;
			close OUTFILE;
		}
	}
}
#write 3_gen_crossover file for maternal
foreach $offspr (keys %offspr_hash){
	if(exists $offspr_mat_hash{$offspr}){
		$numkids = scalar(@{$offspr_mat_hash{$offspr}});
		#for each offspring of a mother who is also a kid in another trio, open the mothers trio file and the offsprings trio file
		@toks=split('_',@{$offspr_hash{$offspr}}[0]);
		$gm_mat = @toks[0];
		$gp_mat = @toks[1];
		foreach $offspring (@{$offspr_mat_hash{$offspr}}){
			$fname_out = "$offspring" . "_crossovers_3gen.csv";
			if(-e $fname_out){
				open (OUTFILE, ">>$fname_out") || die;
			}else{
				open (OUTFILE, ">$fname_out") || die;
			}
			$offspring_trio_filename = "@{$offspr_hash{$offspring}}[0]"."_geno.csv";
			$maternal_trio_filename = "@{$offspr_hash{$offspr}}[0]"."_geno.csv";
			open (INFILE_M, "$maternal_trio_filename") || die;
			open (INFILE_O, "$offspring_trio_filename") || die;
			#read in the headers for each trio file
			$headerline1=<INFILE_M>;
			$headerline2=<INFILE_O>;
			while($linem=<INFILE_M>){
				$lineo=<INFILE_O>;
				chomp $linem;
				chomp $lineo;
				@toksm = split(',',$linem);
				@tokso = split(',',$lineo);
				if((join(',',@toksm[3..9]) =~ /[emu]/) || (join(',',@tokso[3..9]) =~ /[emu]/)){
					next;
				}
				if(@toksm[5] == 1){
					if(@toksm[8] == @tokso[6]){
						print OUTFILE "$gp_mat"."_"."$offspr,";
						print OUTFILE "@tokso[2],";
						$currg = $gp_mat;
						$prev_loc = @tokso[2];
						last;
					}elsif(@toksm[6] == @tokso[6]){
						print OUTFILE "$gm_mat" . "_" ."$offspr,";
						print OUTFILE "@tokso[2],";
						$currg = $gm_mat;
						$prev_loc = @tokso[2];
						last;
					}
				}
			}#found first informative SNP and exited while loop
			$switch_now = 0;
			while($linem=<INFILE_M>){
				$lineo=<INFILE_O>;
				chomp $linem;
				chomp $lineo;
				@toksm = split(',',$linem);
				@tokso = split(',',$lineo);
				if((join(',',@toksm[3..9]) =~ /[emu]/) || (join(',',@tokso[3..9]) =~ /[emu]/)){
					next;
				}
				if(@toksm[5] == 1){
					if(@toksm[6] == @tokso[6]){
						if($currg eq $gm_mat){
							$prev_loc = @tokso[2];
							$switch_now = 0;
							next;
						}else{
							if($switch_now == 0){
								$end_loc = $prev_loc;
								$start_loc = @tokso[2];
							}
							$switch_now++;
							$prev_loc = @tokso[2];
							if($switch_now==@_[0]){
								print OUTFILE "$end_loc\n$gm_mat"."_"."$offspr,$start_loc,";
								$currg = $gm_mat;
								$switch_now=0;
							}
						}
					}elsif(@toksm[8] == @tokso[6]){
						if($currg eq $gp_mat){
							$prev_loc=@tokso[2];
							$switch_now = 0;
							next;
						}else{
							if($switch_now == 0){
								$end_loc=$prev_loc;
								$start_loc=@tokso[2];
							}
							$switch_now++;
							$prev_loc = @tokso[2];
							if($switch_now == @_[0]){
								print OUTFILE "$end_loc\n$gp_mat_";
								print OUTFILE "$offspr,$start_loc,";
								$currg = $gp_mat;
								$switch_now = 0;
							}
						}
					}
				}
			}
			print OUTFILE "$prev_loc\n";
			close INFILE_O;
			close INFILE_M;
			close OUTFILE;
		}
	}
}
}

sub Step3 {
### program phases both parents into two separate files named "[parentid]_phased.csv" and
### writes a file for each offspring of the ranges between each crossover event named "[offsprid]_crossover.csv"

#ARGV[0] is marker requirement for switches
if(@_[0] == 0){
	print STDOUT "A minimum of 1 marker is required to identify a crossover.";
	die;
}
if(@_[0] > 10){
	print STDOUT "You are using >10 markers to identify a crossover. If you observe strange results try lowering the number of markers.";
}

@files = glob("*crossover.csv"); 
unlink @files;
@files2 = glob("*phased.csv"); 
unlink @files2;


###Read in cleared trios file and store hashes with mat and pat as keys
%offspr_hash = ();
%offspr_mat_hash = ();
%offspr_pat_hash = ();

open (INFILE_CLEARED, "Cleared_trios.csv") || die;
while($trioline=<INFILE_CLEARED>){	
	chomp $trioline;
	@toks1 = split('_', $trioline);
	$maternal = @toks1[0];
	$paternal = @toks1[1];
	$offspr = @toks1[2];
	push @mats, $maternal;
	push @pats, $paternal;
	push @kids, $offspr;
	@{$offspr_hash{$offspr}}[0] = $trioline;
	push @{$offspr_pat_hash{$paternal}}, $offspr;
	push @{$offspr_mat_hash{$maternal}}, $offspr;
}

#for each mother, apply minimum recombination model
foreach $currentmom (keys %offspr_mat_hash){
	$numberkids=scalar(@{$offspr_mat_hash{$currentmom}});
	if($numberkids < 3){
		next;
	}
	
	#open outfile for phased mom file
	$mfname = "$currentmom" . "_phased.csv";
	open (M_OUTFILE, ">$mfname") || die;
	
	#open infile and outfiles for kids
	@offspringinfilearray = ();
	$i = 0;
	foreach $offspring (@{$offspr_mat_hash{$currentmom}}){
		$kidinfilename = "@{$offspr_hash{$offspring}}[0]"."_geno.csv";
		$kidoutfilename = "$offspring"."_crossover.csv";
		open ($offspringinfilearray[$i], "$kidinfilename") || die;
		open ($offspringoutfilearray[$i], ">$kidoutfilename") || die;
		$i++;
	}
	print M_OUTFILE "SNP,Chr,Position,$currentmom" . "_1," . "$currentmom" . "_2\n";

	#traverse to first heterozygous maternal SNP
	#read in header line 
	for($i = 0; $i < $numberkids; $i++){	
		$lines[$i] = readline($offspringinfilearray[$i]);
	}
	
	#boolean flag because homozygous maternal SNPs provide no information, need to find first heterozygous SNP
	$start = 0;
	while ($start == 0){
		$minct = 0;
		for($i = 0; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chomp $lines[$i];		
		}
		@toks = split(',', $lines[0]);

		if($toks[3] != 1){
			if($toks[3] =~ /[emu]/){
				$toks[3] = "-";
			}
			elsif($toks[3] == 0){
				$toks[3] = "0";
			}
			elsif($toks[3] == 2){
				$toks[3] = "1";
			}
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[3],$toks[3]\n";
			next;
		}
		else{
			for($i = 0; $i < $numberkids; $i++){
				@toks = split(',', $lines[$i]);
				if ($toks[6] !~ /[emu]/){			
					$minct++;
				}
			}
			if($minct != $numberkids){
				print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";
				next;
			}
			else{
				$start = 1;
			}
		}		
	}
	#initial heterozygous case to determine which chr the offspr are on, we already have the current @lines
	#keep @chr_ary for each offspr's current chr1 or chr2 (choose the chr with the value '1' for each offspr)
	#must have three offspr genotype minimum or skip 
	for($i = 0; $i < $numberkids; $i++){
		@toks = split(',', $lines[$i]);
		if($toks[6] == 1){
			#set to chr1 (index 6 in _geno file), else to chr2 (index 7)
			$chr_ary[$i] = 6;
		}
		else{
			$chr_ary[$i] = 7;
		}
	}
	print M_OUTFILE "$toks[0],$toks[1],$toks[2],1,0\n";	
	#write out the parent chr_1 and beginning location of the first segment for each offspring file
	for($i = 0; $i < $numberkids; $i++){
		$prev_hetero_loc[$i] = $toks[2];
		$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];
		if($chr_ary[$i] == 6){
			print {$offspringoutfilearray[$i]} "$currentmom" . "_1," ."$toks[2],";				
		}
		else{
			print {$offspringoutfilearray[$i]} "$currentmom" . "_2," ."$toks[2],";				
		}
	}
	#begin phasing after initializing parent and offspr's @chr_ary for each offspr's current chr1 or chr2
	while(!(eof $offspringinfilearray[0])){
		$minct = 0;
		#if parent is homozygous, we already know the output
		$lines[0] = readline($offspringinfilearray[0]);
		chop $lines[0];
		@toks = split(',', $lines[0]);	
		if(($toks[3] == 0) || ($toks[3] == 2) || ($toks[3] =~ /[emu]/)){
			if($toks[3] =~ /[emu]/){
				$toks[3] = "-";
			}
			elsif($toks[3] == 0){
				$toks[3] = "0";
			}
			elsif($toks[3] == 2){
				$toks[3] = "1";
			}

			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[3],$toks[3]\n";
			#advance to next line for rest of offspr
			for($i = 1; $i < $numberkids; $i++){	
				$lines[$i] = readline($offspringinfilearray[$i]);
			}
			next;			
		}
		#parent is heterozygous
		#already read in line from first offspr so start at $i = 1, take care of first offspr first
		if ($toks[6] !~ /[emu]/){							
			$val_ary[0] = $toks[$chr_ary[0]];
			$minct++;				
		}
		else{
			$val_ary[0] = "";
		}
		for($i = 1; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chop $lines[$i];
			@toks = split(',', $lines[$i]);
			if ($toks[6] !~ /[emu]/){							
				$val_ary[$i] = $toks[$chr_ary[$i]];
				$minct++;
			}
			else{
				$val_ary[$i] = "";
			}
		}
		if($minct >2){
			#use a popular vote of the "sum" of the 0's or 1's to determine the genotype
			#and to determine if there is a crossover event 
			$sum = 0;
			for($i = 0; $i < $numberkids; $i++){
				$sum = $sum + $val_ary[$i];
				if($val_ary[$i] !~ /[emu]/){
					$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];				
					$prev_hetero_loc[$i] = $toks[2];
				}
			}
			#no crossover event, all 0's or all 1's
			if(($sum == $minct) || ($sum == 0)){
				for($i = 0; $i < $numberkids; $i++){
					@counts[$i]=();
				}			
				if($sum == $minct){
					$m_chr1 = 1;
					$m_chr2 = 0;
				}
				else{
					$m_chr1 = 0;			
					$m_chr2 = 1;
				}						
			}
			#crossover event
			else{
				#count the number of 0's
				$ct = 0;
				for($i = 0; $i < $numberkids; $i++){
					if($val_ary[$i] eq "0"){
						$ct++;
					}
				}
				#count the number of 1's
				$ct2 = 0;
				for($j = 0; $j < $numberkids; $j++){
					if($val_ary[$j] eq "1"){
						$ct2++;
					}
				}
				#decide which allele is the current allele
				if($ct>$ct2){
					#the common allele is a 0
					$m_chr1 = 0;
					$m_chr2 = 1;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "1"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							if($chr_ary[$i] == 6){
								$chr_ary[$i] = 7;
							}
							elsif($chr_ary[$i] == 7){
								$chr_ary[$i] = 6;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 6){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
				elsif($ct2>$ct){
					#the common allele is a 1
					$m_chr1 = 1;
					$m_chr2 = 0;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							@counts[$i] = ();
							#change crossover offspr strand	
							if($chr_ary[$i] == 6){
								$chr_ary[$i] = 7;
							}
							elsif($chr_ary[$i] == 7){
								$chr_ary[$i] = 6;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 6){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}		
						}
					}
				}
				elsif($ct2==$ct){
					#throw out SNP, its ambiguous, what if this happens?
					$m_chr1 = "-";
					$m_chr2 = "-";
					for($i = 0; $i < $numberkids; $i++){
						#print STDOUT "$val_ary[$i]\t";
						#print STDOUT "\n";
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							if($chr_ary[$i] == 6){
								$chr_ary[$i] = 7;
							}
							elsif($chr_ary[$i] == 7){
								$chr_ary[$i] = 6;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 6){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
			}
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$m_chr1,$m_chr2\n";	
		}
		#there weren't at least three offspr with genotype data
		else{
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";	
		}
	}
	#print end of last segment for each offspring (last heterozygous SNP and close the files
	for($i = 0; $i < $numberkids; $i++){
		print {$offspringoutfilearray[$i]} "$prev_hetero_loc[$i]\n";
		close $offspringoutfilearray[$i];
	}
	for($i = 0; $i < $numberkids; $i++){
		close $offspringinfilearray[$i];
	}
	close M_OUTFILE;
}

#for each father, apply minimum recombination model
foreach $currentdad (keys %offspr_pat_hash){
	$numberkids=scalar(@{$offspr_pat_hash{$currentdad}});
	if($numberkids < 3){
		next;
	}
	
	#open outfile for phased dad file
	$mfname = "$currentdad" . "_phased.csv";
	open (M_OUTFILE, ">$mfname") || die;
	
	#open infile and outfiles for kids
	@offspringinfilearray = ();
	$i = 0;
	foreach $offspring (@{$offspr_pat_hash{$currentdad}}){
		$kidinfilename = "@{$offspr_hash{$offspring}}[0]"."_geno.csv";
		$kidoutfilename = "$offspring"."_crossover.csv";
		open ($offspringinfilearray[$i], "$kidinfilename") || die;
		open ($offspringoutfilearray[$i], ">>$kidoutfilename") || die;
		$i++;
	}
	print M_OUTFILE "SNP,Chr,Position,$currentdad" . "_1," . "$currentdad" . "_2\n";

	#traverse to first heterozygous maternal SNP
	#read in header line 
	for($i = 0; $i < $numberkids; $i++){	
		$lines[$i] = readline($offspringinfilearray[$i]);
	}
	
	#boolean flag because homozygous maternal SNPs provide no information, need to find first heterozygous SNP
	$start = 0;
	while ($start == 0){
		$minct = 0;
		for($i = 0; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chomp $lines[$i];		
		}
		@toks = split(',', $lines[0]);
		if($toks[4] != 1){
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";
			next;
		}
		else{
			for($i = 0; $i < $numberkids; $i++){
				@toks = split(',', $lines[$i]);
				if ($toks[8] !~ /[emu]/){			
					$minct++;
				}
			}
			if($minct != $numberkids){
				if($toks[4] =~ /[emu]/){
					$toks[4] = "-";
				}
				elsif($toks[4] == 0){
					$toks[4] = "0";
				}
				elsif($toks[4] == 2){
					$toks[4] = "1";
				}
				print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[4],$toks[4]\n";
				next;
			}
			else{
				$start = 1;
			}
		}		
	}
	#initial heterozygous case to determine which chr the offspr are on, we already have the current @lines
	#keep @chr_ary for each offspr's current chr1 or chr2 (choose the chr with the value '1' for each offspr)
	#must have three offspr genotype minimum or skip 
	for($i = 0; $i < $numberkids; $i++){
		@toks = split(',', $lines[$i]);
		if($toks[8] == 1){
			#set to chr1 (index $mind1 in _geno file), else to chr2 (index $mind2)
			$chr_ary[$i] = 8;
		}
		else{
			$chr_ary[$i] = 9;
		}
	}
	print M_OUTFILE "$toks[0],$toks[1],$toks[2],1,0\n";	
	#write out the parent chr_1 and beginning location of the first segment for each offspring file
	for($i = 0; $i < $numberkids; $i++){
		$prev_hetero_loc[$i] = $toks[2];
		$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];
		if($chr_ary[$i] == 8){
			print {$offspringoutfilearray[$i]} "$currentdad" . "_1," ."$toks[2],";				
		}
		else{
			print {$offspringoutfilearray[$i]} "$currentdad" . "_2," ."$toks[2],";				
		}
	}
	#begin phasing after initializing parent and offspr's @chr_ary for each offspr's current chr1 or chr2
	while(!(eof $offspringinfilearray[0])){
		$minct = 0;
		#if parent is homozygous, we already know the output
		$lines[0] = readline($offspringinfilearray[0]);
		chop $lines[0];
		@toks = split(',', $lines[0]);	
		if(($toks[4] == 0) || ($toks[4] == 2) || ($toks[4] =~ /[emu]/)){
			if($toks[4] =~ /[emu]/){
				$toks[4] = "-";
			}
			elsif($toks[4] == 0){
				$toks[4] = "0";
			}
			elsif($toks[4] == 2){
				$toks[4] = "1";
			}

			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[4],$toks[4]\n";
			#advance to next line for rest of offspr
			for($i = 1; $i < $numberkids; $i++){	
				$lines[$i] = readline($offspringinfilearray[$i]);
			}
			next;			
		}
		#parent is heterozygous
		#already read in line from first offspr so start at $i = 1, take care of first offspr first
		if ($toks[8] !~ /[emu]/){							
			$val_ary[0] = $toks[$chr_ary[0]];
			$minct++;				
		}
		else{
			$val_ary[0] = "";
		}
		for($i = 1; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chop $lines[$i];
			@toks = split(',', $lines[$i]);
			if ($toks[8] !~ /[emu]/){							
				$val_ary[$i] = $toks[$chr_ary[$i]];
				$minct++;
			}
			else{
				$val_ary[$i] = "";
			}
		}
		if($minct >2){
			#use a popular vote of the "sum" of the 0's or 1's to determine the genotype
			#and to determine if there is a crossover event 
			$sum = 0;
			for($i = 0; $i < $numberkids; $i++){
				$sum = $sum + $val_ary[$i];
				if($val_ary[$i] !~ /[emu]/){
					$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];				
					$prev_hetero_loc[$i] = $toks[2];
				}
			}
			#no crossover event, all 0's or all 1's
			if(($sum == $minct) || ($sum == 0)){
				for($i = 0; $i < $numberkids; $i++){
					@counts[$i]=();
				}			
				if($sum == $minct){
					$m_chr1 = 1;
					$m_chr2 = 0;
				}
				else{
					$m_chr1 = 0;			
					$m_chr2 = 1;
				}						
			}
			#crossover event
			else{
				#count the number of 0's
				$ct = 0;
				for($i = 0; $i < $numberkids; $i++){
					if($val_ary[$i] eq "0"){
						$ct++;
					}
				}
				#count the number of 1's
				$ct2 = 0;
				for($j = 0; $j < $numberkids; $j++){
					if($val_ary[$j] eq "1"){
						$ct2++;
					}
				}
				#decide which allele is the current allele
				if($ct>$ct2){
					#the common allele is a 0
					$m_chr1 = 0;
					$m_chr2 = 1;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "1"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							if($chr_ary[$i] == 8){
								$chr_ary[$i] = 9;
							}
							elsif($chr_ary[$i] == 9){
								$chr_ary[$i] = 8;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 8){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
				elsif($ct2>$ct){
					#the common allele is a 1
					$m_chr1 = 1;
					$m_chr2 = 0;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							@counts[$i] = ();
							#change crossover offspr strand	
							if($chr_ary[$i] == 8){
								$chr_ary[$i] = 9;
							}
							elsif($chr_ary[$i] == 9){
								$chr_ary[$i] = 8;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 8){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}		
						}
					}
				}
				elsif($ct2==$ct){
					#throw out SNP, its ambiguous, what if this happens?
					$m_chr1 = "-";
					$m_chr2 = "-";
					for($i = 0; $i < $numberkids; $i++){
						#print STDOUT "$val_ary[$i]\t";
						#print STDOUT "\n";
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						#print {$offspringoutfilearray[$i]} "Entered common allele parsimonious = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray[$i]} "Entered parsimonious switch count 1\n";
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							#print {$offspringoutfilearray[$i]} "Entered parsimony switch\n";
							if($chr_ary[$i] == 8){
								$chr_ary[$i] = 9;
							}
							elsif($chr_ary[$i] == 9){
								$chr_ary[$i] = 8;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 8){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
			}
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$m_chr1,$m_chr2\n";	
		}
		#there weren't at least three offspr with genotype data
		else{
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";	
		}
	}
	#print end of last segment for each offspring (last heterozygous SNP and close the files
	for($i = 0; $i < $numberkids; $i++){
		print {$offspringoutfilearray[$i]} "$prev_hetero_loc[$i]\n";
		close $offspringoutfilearray[$i];
	}
	for($i = 0; $i < $numberkids; $i++){
		close $offspringinfilearray[$i];
	}
	close M_OUTFILE;
}
}
sub Step4 {
### program phases both parents into two separate files named "[parentid]_phased.csv" and
### writes a file for each offspring of the ranges between each crossover event named "[offsprid]_crossover.csv"

#ARGV[0] is marker requirement for switches
if(@_[0] == 0){
	print STDOUT "A minimum of 1 marker is required to identify a crossover.";
	die;
}
if(@_[0] > 10){
	print STDOUT "You are using >10 markers to identify a crossover. If you observe strange results try lowering the number of markers.";
}

###Read in cleared trios file and store hashes with mat and pat as keys
%offspr_hash = ();
%offspr_mat_hash = ();
%offspr_pat_hash = ();

open (INFILE_CLEARED, "Cleared_trios.csv") || die;
while($trioline=<INFILE_CLEARED>){	
	chomp $trioline;
	@toks1 = split('_', $trioline);
	$maternal = @toks1[0];
	$paternal = @toks1[1];
	$offspr = @toks1[2];
	push @mats, $maternal;
	push @pats, $paternal;
	push @kids, $offspr;
	@{$offspr_hash{$offspr}}[0] = $trioline;
	push @{$offspr_pat_hash{$paternal}}, $offspr;
	push @{$offspr_mat_hash{$maternal}}, $offspr;
}

#for each mother, apply minimum recombination model
foreach $currentmom (keys %offspr_mat_hash){
	$numberkids=scalar(@{$offspr_mat_hash{$currentmom}});
	if($numberkids != 2){
		next;
	}
	
	#open outfile for phased mom file
	$mfname = "$currentmom" . "_phased.csv";
	open (M_OUTFILE, ">$mfname") || die;
	
	#open infile and outfiles for kids
	@offspringinfilearray = ();
	$i = 0;
	foreach $offspring (@{$offspr_mat_hash{$currentmom}}){
		$kidinfilename = "@{$offspr_hash{$offspring}}[0]"."_geno.csv";
		$kidoutfilename = "$offspring"."_crossover.csv";
		open ($offspringinfilearray[$i], "$kidinfilename") || die;
		open ($offspringoutfilearray[$i], ">$kidoutfilename") || die;
		$i++;
	}
	print M_OUTFILE "SNP,Chr,Position,$currentmom" . "_1," . "$currentmom" . "_2\n";

	#traverse to first heterozygous maternal SNP
	#read in header line 
	for($i = 0; $i < $numberkids; $i++){	
		$lines[$i] = readline($offspringinfilearray[$i]);
	}
	
	#boolean flag because homozygous maternal SNPs provide no information, need to find first heterozygous SNP
	$start = 0;
	while ($start == 0){
		$minct = 0;
		for($i = 0; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chomp $lines[$i];		
		}
		@toks = split(',', $lines[0]);

		if($toks[3] != 1){
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";
			next;
		}
		else{
			for($i = 0; $i < $numberkids; $i++){
				@toks = split(',', $lines[$i]);
				if ($toks[6] !~ /[emu]/){			
					$minct++;
				}
			}
			if($minct != $numberkids){
				if($toks[3] =~ /[emu]/){
					$toks[3] = "-";
				}
				elsif($toks[3] == 0){
					$toks[3] = "0";
				}
				elsif($toks[3] == 2){
					$toks[3] = "1";
				}
				print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[3],$toks[3]\n";
				next;
			}
			else{
				$start = 1;
			}
		}		
	}
	#initial heterozygous case to determine which chr the offspr are on, we already have the current @lines
	#keep @chr_ary for each offspr's current chr1 or chr2 (choose the chr with the value '1' for each offspr)
	#must have three offspr genotype minimum or skip 
	for($i = 0; $i < $numberkids; $i++){
		@toks = split(',', $lines[$i]);
		if($toks[6] == 1){
			#set to chr1 (index 6 in _geno file), else to chr2 (index 7)
			$chr_ary[$i] = 6;
		}
		else{
			$chr_ary[$i] = 7;
		}
	}
	print M_OUTFILE "$toks[0],$toks[1],$toks[2],1,0\n";	
	#write out the parent chr_1 and beginning location of the first segment for each offspring file
	for($i = 0; $i < $numberkids; $i++){
		$prev_hetero_loc[$i] = $toks[2];
		$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];
		if($chr_ary[$i] == 6){
			print {$offspringoutfilearray[$i]} "$currentmom" . "_1," ."$toks[2],";				
		}
		else{
			print {$offspringoutfilearray[$i]} "$currentmom" . "_2," ."$toks[2],";				
		}
	}
	#begin phasing after initializing parent and offspr's @chr_ary for each offspr's current chr1 or chr2
	while(!(eof $offspringinfilearray[0])){
		$minct = 0;
		#if parent is homozygous, we already know the output
		$lines[0] = readline($offspringinfilearray[0]);
		chop $lines[0];
		@toks = split(',', $lines[0]);	
		if(($toks[3] == 0) || ($toks[3] == 2) || ($toks[3] =~ /[emu]/)){
			if($toks[3] =~ /[emu]/){
				$toks[3] = "-";
			}
			elsif($toks[3] == 0){
				$toks[3] = "0";
			}
			elsif($toks[3] == 2){
				$toks[3] = "1";
			}

			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[3],$toks[3]\n";
			#advance to next line for rest of offspr
			for($i = 1; $i < $numberkids; $i++){	
				$lines[$i] = readline($offspringinfilearray[$i]);
			}
			next;			
		}
		#parent is heterozygous
		#already read in line from first offspr so start at $i = 1, take care of first offspr first
		if ($toks[6] !~ /[emu]/){							
			$val_ary[0] = $toks[$chr_ary[0]];
			$minct++;				
		}
		else{
			$val_ary[0] = "";
		}
		for($i = 1; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chop $lines[$i];
			@toks = split(',', $lines[$i]);
			if ($toks[6] !~ /[emu]/){							
				$val_ary[$i] = $toks[$chr_ary[$i]];
				$minct++;
			}
			else{
				$val_ary[$i] = "";
			}
		}
		if($minct == 2){
			#use a popular vote of the "sum" of the 0's or 1's to determine the genotype
			#and to determine if there is a crossover event 
			$sum = 0;
			for($i = 0; $i < $numberkids; $i++){
				$sum = $sum + $val_ary[$i];
				if($val_ary[$i] !~ /[emu]/){
					$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];				
					$prev_hetero_loc[$i] = $toks[2];
				}
			}
			#no crossover event, all 0's or all 1's
			if(($sum == $minct) || ($sum == 0)){
				for($i = 0; $i < $numberkids; $i++){
					@counts[$i]=0;
				}			
				if($sum == $minct){
					$m_chr1 = 1;
					$m_chr2 = 0;
				}
				else{
					$m_chr1 = 0;			
					$m_chr2 = 1;
				}						
			}
			#crossover event
			else{
				#count the number of 0's
				$ct = 0;
				for($i = 0; $i < $numberkids; $i++){
					if($val_ary[$i] eq "0"){
						$ct++;
					}
				}
				#count the number of 1's
				$ct2 = 0;
				for($j = 0; $j < $numberkids; $j++){
					if($val_ary[$j] eq "1"){
						$ct2++;
					}
				}
				#decide which allele is the current allele
				if($ct>$ct2){
					#the common allele is a 0
					$m_chr1 = 0;
					$m_chr2 = 1;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "1"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							if($chr_ary[$i] == 6){
								$chr_ary[$i] = 7;
							}
							elsif($chr_ary[$i] == 7){
								$chr_ary[$i] = 6;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 6){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
				elsif($ct2>$ct){
					#the common allele is a 1
					$m_chr1 = 1;
					$m_chr2 = 0;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							@counts[$i] = ();
							#change crossover offspr strand	
							if($chr_ary[$i] == 6){
								$chr_ary[$i] = 7;
							}
							elsif($chr_ary[$i] == 7){
								$chr_ary[$i] = 6;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 6){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}		
						}
					}
				}
				elsif($ct2==$ct){
					#throw out SNP, its ambiguous, what if this happens?
					$m_chr1 = "-";
					$m_chr2 = "-";
					for($i = 0; $i < $numberkids; $i++){
						#print STDOUT "$val_ary[$i]\t";
						#print STDOUT "\n";
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							if($chr_ary[$i] == 6){
								$chr_ary[$i] = 7;
							}
							elsif($chr_ary[$i] == 7){
								$chr_ary[$i] = 6;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 6){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
			}
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$m_chr1,$m_chr2\n";	
		}
		#there weren't at least three offspr with genotype data
		else{
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";	
		}
	}
	#print end of last segment for each offspring (last heterozygous SNP and close the files
	for($i = 0; $i < $numberkids; $i++){
		print {$offspringoutfilearray[$i]} "$prev_hetero_loc[$i]\n";
		close $offspringoutfilearray[$i];
	}
	for($i = 0; $i < $numberkids; $i++){
		close $offspringinfilearray[$i];
	}
	close M_OUTFILE;
}

#for each father, apply minimum recombination model
foreach $currentdad (keys %offspr_pat_hash){
	$numberkids=scalar(@{$offspr_pat_hash{$currentdad}});
	if($numberkids != 2){
		next;
	}
	
	#open outfile for phased dad file
	$mfname = "$currentdad" . "_phased.csv";
	open (M_OUTFILE, ">$mfname") || die;
	
	#open infile and outfiles for kids
	@offspringinfilearray = ();
	$i = 0;
	foreach $offspring (@{$offspr_pat_hash{$currentdad}}){
		$kidinfilename = "@{$offspr_hash{$offspring}}[0]"."_geno.csv";
		$kidoutfilename = "$offspring"."_crossover.csv";
		open ($offspringinfilearray[$i], "$kidinfilename") || die;
		open ($offspringoutfilearray[$i], ">>$kidoutfilename") || die;
		$i++;
	}
	print M_OUTFILE "SNP,Chr,Position,$currentdad" . "_1," . "$currentdad" . "_2\n";

	#traverse to first heterozygous maternal SNP
	#read in header line 
	for($i = 0; $i < $numberkids; $i++){	
		$lines[$i] = readline($offspringinfilearray[$i]);
	}
	
	#boolean flag because homozygous maternal SNPs provide no information, need to find first heterozygous SNP
	$start = 0;
	while ($start == 0){
		$minct = 0;
		for($i = 0; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chomp $lines[$i];		
		}
		@toks = split(',', $lines[0]);
		if($toks[4] != 1){
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";
			next;
		}
		else{
			for($i = 0; $i < $numberkids; $i++){
				@toks = split(',', $lines[$i]);
				if ($toks[8] !~ /[emu]/){			
					$minct++;
				}
			}
			if($minct != $numberkids){
				if($toks[4] =~ /[emu]/){
					$toks[4] = "-";
				}
				elsif($toks[4] == 0){
					$toks[4] = "0";
				}
				elsif($toks[4] == 2){
					$toks[4] = "1";
				}
				print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[4],$toks[4]\n";
				next;
			}
			else{
				$start = 1;
			}
		}		
	}
	#initial heterozygous case to determine which chr the offspr are on, we already have the current @lines
	#keep @chr_ary for each offspr's current chr1 or chr2 (choose the chr with the value '1' for each offspr)
	#must have three offspr genotype minimum or skip 
	for($i = 0; $i < $numberkids; $i++){
		@toks = split(',', $lines[$i]);
		if($toks[8] == 1){
			#set to chr1 (index $mind1 in _geno file), else to chr2 (index $mind2)
			$chr_ary[$i] = 8;
		}
		else{
			$chr_ary[$i] = 9;
		}
	}
	print M_OUTFILE "$toks[0],$toks[1],$toks[2],1,0\n";	
	#write out the parent chr_1 and beginning location of the first segment for each offspring file
	for($i = 0; $i < $numberkids; $i++){
		$prev_hetero_loc[$i] = $toks[2];
		$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];
		if($chr_ary[$i] == 8){
			print {$offspringoutfilearray[$i]} "$currentdad" . "_1," ."$toks[2],";				
		}
		else{
			print {$offspringoutfilearray[$i]} "$currentdad" . "_2," ."$toks[2],";				
		}
	}
	#begin phasing after initializing parent and offspr's @chr_ary for each offspr's current chr1 or chr2
	while(!(eof $offspringinfilearray[0])){
		$minct = 0;
		#if parent is homozygous, we already know the output
		$lines[0] = readline($offspringinfilearray[0]);
		chop $lines[0];
		@toks = split(',', $lines[0]);	
		if(($toks[4] == 0) || ($toks[4] == 2) || ($toks[4] =~ /[emu]/)){
			if($toks[4] =~ /[emu]/){
				$toks[4] = "-";
			}
			elsif($toks[4] == 0){
				$toks[4] = "0";
			}
			elsif($toks[4] == 2){
				$toks[4] = "1";
			}

			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$toks[4],$toks[4]\n";
			#advance to next line for rest of offspr
			for($i = 1; $i < $numberkids; $i++){	
				$lines[$i] = readline($offspringinfilearray[$i]);
			}
			next;			
		}
		#parent is heterozygous
		#already read in line from first offspr so start at $i = 1, take care of first offspr first
		if ($toks[8] !~ /[emu]/){							
			$val_ary[0] = $toks[$chr_ary[0]];
			$minct++;				
		}
		else{
			$val_ary[0] = "";
		}
		for($i = 1; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chop $lines[$i];
			@toks = split(',', $lines[$i]);
			if ($toks[8] !~ /[emu]/){							
				$val_ary[$i] = $toks[$chr_ary[$i]];
				$minct++;
			}
			else{
				$val_ary[$i] = "";
			}
		}
		if($minct == 2){
			#use a popular vote of the "sum" of the 0's or 1's to determine the genotype
			#and to determine if there is a crossover event 
			$sum = 0;
			for($i = 0; $i < $numberkids; $i++){
				$sum = $sum + $val_ary[$i];
				if($val_ary[$i] !~ /[emu]/){
					$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];				
					$prev_hetero_loc[$i] = $toks[2];
				}
			}
			#no crossover event, all 0's or all 1's
			if(($sum == $minct) || ($sum == 0)){
				for($i = 0; $i < $numberkids; $i++){
					@counts[$i]=();
				}			
				if($sum == $minct){
					$m_chr1 = 1;
					$m_chr2 = 0;
				}
				else{
					$m_chr1 = 0;			
					$m_chr2 = 1;
				}						
			}
			#crossover event
			else{
				#count the number of 0's
				$ct = 0;
				for($i = 0; $i < $numberkids; $i++){
					if($val_ary[$i] eq "0"){
						$ct++;
					}
				}
				#count the number of 1's
				$ct2 = 0;
				for($j = 0; $j < $numberkids; $j++){
					if($val_ary[$j] eq "1"){
						$ct2++;
					}
				}
				#decide which allele is the current allele
				if($ct>$ct2){
					#the common allele is a 0
					$m_chr1 = 0;
					$m_chr2 = 1;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "1"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							if($chr_ary[$i] == 8){
								$chr_ary[$i] = 9;
							}
							elsif($chr_ary[$i] == 9){
								$chr_ary[$i] = 8;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 8){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
				elsif($ct2>$ct){
					#the common allele is a 1
					$m_chr1 = 1;
					$m_chr2 = 0;
					for($i = 0; $i < $numberkids; $i++){
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
						}
						if(@_[0]==$currcount1){
							@counts[$i] = ();
							#change crossover offspr strand	
							if($chr_ary[$i] == 8){
								$chr_ary[$i] = 9;
							}
							elsif($chr_ary[$i] == 9){
								$chr_ary[$i] = 8;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 8){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}		
						}
					}
				}
				elsif($ct2==$ct){
					#throw out SNP, its ambiguous, what if this happens?
					$m_chr1 = "-";
					$m_chr2 = "-";
					for($i = 0; $i < $numberkids; $i++){
						#print STDOUT "$val_ary[$i]\t";
						#print STDOUT "\n";
						if($val_ary[$i] eq "0"){
							@counts[$i]++;
						}
					}
					for($i = 0; $i < $numberkids; $i++){
						$currcount1 = @counts[$i];
						$currcountone=1;
						#print {$offspringoutfilearray[$i]} "Entered common allele parsimonious = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray[$i]} "Entered parsimonious switch count 1\n";
						}
						if(@_[0]==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							#print {$offspringoutfilearray[$i]} "Entered parsimony switch\n";
							if($chr_ary[$i] == 8){
								$chr_ary[$i] = 9;
							}
							elsif($chr_ary[$i] == 9){
								$chr_ary[$i] = 8;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == 8){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_1,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentdad" . "_2,";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i],";
							}	
						}
					}
				}
			}
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],$m_chr1,$m_chr2\n";	
		}
		#there weren't at least three offspr with genotype data
		else{
			print M_OUTFILE "$toks[0],$toks[1],$toks[2],-,-\n";	
		}
	}
	#print end of last segment for each offspring (last heterozygous SNP and close the files
	for($i = 0; $i < $numberkids; $i++){
		print {$offspringoutfilearray[$i]} "$prev_hetero_loc[$i]\n";
		close $offspringoutfilearray[$i];
	}
	for($i = 0; $i < $numberkids; $i++){
		close $offspringinfilearray[$i];
	}
	close M_OUTFILE;
}

}
sub Step5 {
#open list of all trios that have the offspring phased
open (INFILE_CLEARED, "Cleared_trios.csv") || die;

#open founder outfile
open (OUTFILE, ">Nonfounder_matrix.csv") || die;

@trios_ary = <INFILE_CLEARED>;
$filelen = scalar(@trios_ary);

%founder_hash = ();
%second_gen_hash = ();  
@founder_kids = ();
%parents_kids = ();
@tmp_ary = ();
%offspr_hash = ();

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

#find all children of each founder id's by traversing @trios_ary and storing trio file in @founder_kid
#maternal
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];
	if(exists $offspr_hash{$maternal}){
		#not a founder
		next;
	}
	else{
		#is a founder		
		push(@founder_kids, $offspr);
		$founder_hash{$maternal} = 1;
		@tmp_ary = @{$second_gen_hash{$maternal}};
		push(@tmp_ary, $trioline);
		$second_gen_hash{$maternal} = [@tmp_ary];
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
		#not a founder
		next;
	}
	else{
		#is a founder		
		push(@founder_kids, $offspr);
		$founder_hash{$paternal} = 1;
		@tmp_ary = @{$second_gen_hash{$paternal}};
		push(@tmp_ary, $trioline);
		$second_gen_hash{$paternal} = [@tmp_ary];
	}
}
print OUTFILE "Individ,Founding_chr,Start,End,Path\n";
foreach $founder (keys %second_gen_hash){
	@tmp_ary = @{$second_gen_hash{$founder}};
	$len = @tmp_ary;
	for($i = 0; $i < $len; $i++){
		@triotoks = split('_', $tmp_ary[$i]);
		$maternal = @triotoks[0];
		$paternal = @triotoks[1];
		$offspr = @triotoks[2];
		if($len >2){	
			$fnamecross = "$offspr" . "_crossover.csv";
			if (-e $fnamecross){
				open (INFILE_CROSS, "$fnamecross");
				while($line = <INFILE_CROSS>){
					chop $line;
					@toks = split(",", $line);
					if ($toks[0] =~ /$founder/){
						print OUTFILE "$offspr,";
						print OUTFILE "$toks[0],$toks[1],$toks[2],";
						$path = "$offspr" . "_" . "$toks[0]";
						print OUTFILE "$path\n";
					}
				}
				close INFILE_CROSS;
			}
		}
		elsif($len == 2){
			$fnamecross = "$offspr" . "_crossover.csv";
			if (-e $fnamecross){
				open (INFILE_CROSS, "$fnamecross");
				while($line = <INFILE_CROSS>){
					chop $line;
					@toks = split(',', $line);
					if ($toks[0] =~ /$founder/){
						print OUTFILE "$offspr,";
						print OUTFILE "$toks[0],$toks[1],$toks[2],";
						$path = "$offspr" . "_" . "$toks[0]";
						print OUTFILE "$path\n";
					}
				}
				close INFILE_CROSS;
			}
		}
		else{
			#this founder only has either one or two children since no crossover file exists
			$fname_geno = "$tmp_ary[0]" . "_geno.csv";
			open (INFILE_GENO, $fname_geno);
			$headerline = <INFILE_GENO>;
			@lines = <INFILE_GENO>;
			$num_lines = @lines;						
			print OUTFILE "$offspr,";
			print OUTFILE "$founder" . "_" . "1,";
			@toks = split(",", $lines[0]);
			$start = $toks[2];
			print OUTFILE "$start,";
			@toks = split(",", $lines[$num_lines-1]);
			$end = $toks[2];
			print OUTFILE "$end,";
			print OUTFILE "$offspr" ."_". "$founder" . "_1\n";
			close INFILE_GENO;			
		}
	}	
}
close OUTFILE;

########store all children in parents_kids hash, with key as mother or father (unless they are a founder), and value is an array of the children trios
%parents_kids = ();
for($i = 0; $i < $filelen; $i++){	
	#maternal
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];
	if(exists $founder_hash{$maternal}){
		#is a founder
		next;
	}
	else{
		@tmp_ary = @{$parents_kids{$maternal}};
		push(@tmp_ary, $trioline);
		$parents_kids{$maternal} = [@tmp_ary];	
	}
}

for($i = 0; $i < $filelen; $i++){	
	#paternal
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];
	if(exists $founder_hash{$paternal}){
		#is a founder
		next;
	}
	else{
		@tmp_ary = @{$parents_kids{$paternal}};
		push(@tmp_ary, $trioline);
		$parents_kids{$paternal} = [@tmp_ary];	
	}
}

foreach $kid (@founder_kids){
	&print_sub($kid);
}
}

sub print_sub{
	my $this_kid = $_[0];
	if(exists $parents_kids{$this_kid}){
		my @tmp_ary = @{$parents_kids{$this_kid}};
		my $len = @tmp_ary;
		my $i;
		for ($i = 0; $i < $len; $i++){			
			my @triotoks = split('_', $tmp_ary[$i]);
			$maternal = @triotoks[0];
			$paternal = @triotoks[1];
			my $offspr = @triotoks[2];
			print "$offspr\n";
			&chop_sub($offspr);
			&print_sub($offspr);
		}
	}
}

sub chop_sub{
	my $this_offspr = $_[0];
	my $fname = "$this_offspr" . "_crossovers_3gen.csv";
	open (INFILETHREE, $fname) || die;
	open (OUTFILE, ">>Nonfounder_matrix.csv") || die;	
	my @lines3gen = <INFILETHREE>;
	my $line3;
	my $i;
	my $len = @lines3gen;
	for($i = 0; $i < $len; $i++){
		$line3 = $lines3gen[$i];
		chop $line3;
		my @toks = split(',', $line3);
		my $pair3 = $toks[0];
		@toks2 = split('_', $toks[0]);
		#rearranging order of 3gen file ids to match nonfounder matrix
		my $id1 = $toks2[1];
		my $id2 = $toks2[0];
		open (INFILEMATRIX, "Nonfounder_matrix.csv") || die;
		my $linemat = <INFILEMATRIX>;
		while($linemat = <INFILEMATRIX>){
			chop $linemat;
			my @toksmat = split(',', $linemat);
			@toksmat2 = split('_', $toksmat[4]);
			my $id1mat = $toksmat2[0];
			my $id2mat = $toksmat2[1];
			my @tmp_ary;
			my $b3 = $toks[1];
			my $c3 = $toks[2];
			my $d3 = $toksmat[2];
			my $e3 = $toksmat[3];			
			my $newpath = "$this_offspr" . "_" . "$toksmat[4]";
			if(($id1 eq $id1mat) && ($id2 eq $id2mat)){
				if($b3 > $e3){
					next;
				}
				elsif($c3 < $d3){
					next;
				}
				elsif($b3 < $d3){
					if($c3 == $d3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$d3,$d3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 < $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$d3,$c3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 == $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$d3,$c3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 > $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$d3,$e3,";
						print OUTFILE "$newpath\n";
					}
					else{
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "Error,Error,";
						print OUTFILE "$newpath\n";
					}
				}
				elsif($b3 == $d3){
					if($c3 < $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$b3,$c3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 == $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$b3,$c3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 > $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$b3,$e3,";
						print OUTFILE "$newpath\n";
					}
					else{
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "Error,Error,";
						print OUTFILE "$newpath\n";
					}
				}
				elsif($b3 > $d3){
					if($c3 < $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$b3,$c3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 == $e3){
						print OUTFILE "$this_offspr,$toksmat[1],";
						print OUTFILE "$b3,$c3,";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 > $e3){
						if($b3 == $e3){
							print OUTFILE "$this_offspr,$toksmat[1],";
							print OUTFILE "$b3,$e3,";
							print OUTFILE "$newpath\n";
						}
						else{
							print OUTFILE "$this_offspr,$toksmat[1],";
							print OUTFILE "$b3,$e3,";
							print OUTFILE "$newpath\n";
						}
					}
				}
				else{
					print OUTFILE "$this_offspr,$toksmat[1],";
					print OUTFILE "ERROR,$newpath\n";
				}	
			}			
		}
		close INFILEMATRIX;	
	}
	close INFILETHREE;	
	close OUTFILE;
}

sub Step6 {
my $file = "Nonfounder_matrix.csv";
my %seen = ();
{
   @ARGV = ($file);
   $^I = '.bac';
   while(<>){
      $seen{$_}++;
      next if $seen{$_} > 1;
      print;
   }
}
}

1;