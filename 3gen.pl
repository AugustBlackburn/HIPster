#!/usr/local/bin/perl


#ARGV[0] is pedigree number, $ARGV[1] is the number of markers required for a switch to be considered "real"
$ped = $ARGV[0];
$real_size = $ARGV[1];
print STDOUT "The number of markers required for a switch is $real_size.\n";
if($real_size==0){
	print STDOUT "Please enter a value greater than 0 for the number of markers necessary for switches to be considered real\n";
	die;
}
if($real_size > 10){
	print STDOUT "You have indicated that more than 10 markers are necessary to support a recombination event. Consider using less markers if you obtain odd results.\n";
}

#open list of all trios that have the offspring phased, store in @trios_ary
$filename_cleared = "Cleared_Ped_" . "$ped" . ".txt";
open (INFILE_CLEARED, $filename_cleared) || die;

@trios_ary = <INFILE_CLEARED>;
$filelen = @trios_ary;

###traverse @trios_ary and store hashes with mat and pat as keys, trio file(s) as value(s) if the mat or pat is also an offspring
%offspr_hash = ();
%offspr_mat_hash = ();
%offspr_pat_hash = ();
@tmp_ary = ();

#first collect all offspr into $offspr_hash
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];	
	chop $trioline;
	@triotoks = split('_', $trioline);
	$offspr = @triotoks[2];

	push(@tmp_ary, $trioline);
	$offspr_hash{$offspr} = [@tmp_ary];
	@tmp_ary = ();
}
		
#check all mats if they exist in $offspr_hash
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];
	chop $trioline;
	@triotoks = split('_', $trioline);	
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];

	if(exists $offspr_hash{$maternal}){		
		@offspr_ary = @{$offspr_hash{$maternal}};		
		push(@offspr_ary, $trioline);
	
		@tmp_ary = @{$offspr_mat_hash{$maternal}};
			
		push(@tmp_ary, $offspr_ary[0]);
		push(@tmp_ary, $offspr_ary[1]);
		#add this $trioline to mat_triohash
		$offspr_mat_hash{$maternal} = [@tmp_ary];
	}
}

#check all pats if they exist in $offspr_hash
for($i = 0; $i < $filelen; $i++){	
	$trioline = $trios_ary[$i];
	chop $trioline;
	@triotoks = split('_', $trioline);	
	$maternal = @triotoks[0];
	$paternal = @triotoks[1];
	$offspr = @triotoks[2];

	if(exists $offspr_hash{$paternal}){		
		@offspr_ary = @{$offspr_hash{$paternal}};		
		push(@offspr_ary, $trioline);
	
		@tmp_ary = @{$offspr_pat_hash{$paternal}};		
		push(@tmp_ary, $offspr_ary[0]);
		push(@tmp_ary, $offspr_ary[1]);
		#add this $trioline to pat_triohash
		$offspr_pat_hash{$paternal} = [@tmp_ary];
	}
}
	
@tmp_ary = ();

#write 3_gen_crossover file for paternal
foreach $offspr (keys %offspr_hash){
	if(exists $offspr_pat_hash{$offspr}){
		#write 3_gen_crossover file for paternal part
		@tmp_ary = @{$offspr_pat_hash{$offspr}};

	    $numtrios = @tmp_ary;

	    for($i = 0; $i < $numtrios; $i = $i+2){	
		$trioline = $tmp_ary[$i+1];
		chop $line;
		@triotoks = split('_', $trioline);
		$paternal = $triotoks[1];
		$offspring = $triotoks[2];

		$fname_out = "$offspring" . "_crossovers_3gen.txt";
		open (OUTFILE, ">$fname_out") || die;

		#get paternal's mat and pat (gp_pat, gm_pat)
		$line = $tmp_ary[$i];
		chop $line;
		@toks = split('_', $line);
		$gm_pat = $toks[0];
		$gp_pat = $toks[1];




		$fname_pat = "$tmp_ary[$i]" . "_geno.txt";
		$fname_offspr = "$tmp_ary[$i+1]" . "_geno.txt";
		
		open (INFILE_P, "$fname_pat") || die;
		open (INFILE_O, "$fname_offspr") || die;

		#traverse to first crossover unambiguous chr
		#read in header line from both files
		$line = <INFILE_P>; $line = <INFILE_O>;

		#set index for inherited from gmat or gpat for middle paternal
		$pind1 = 6;
		$pind2 = 8;
		$pindoff = 8;

		$currg = "";

		#boolean flag to first informative SNP
		$start = 0;
		while($start == 0){
			$linep = <INFILE_P>;
			chop $linep;
			@toksp = split('\t', $linep);
			$lineo = <INFILE_O>;
			chop $lineo;
			@tokso = split('\t', $lineo);

			if(($linep =~ /[emu]/) || ($lineo =~ /[emu]/)){				
				next;
			}

			if($toksp[$pind1] == $tokso[$pindoff]){
				if($toksp[$pind2] == $tokso[$pindoff]){				
					next;
				}
				else{
					print OUTFILE "$gm_pat" . "_" ."$paternal\t";
					print OUTFILE "$tokso[2]\t";
					$currg = $gm_pat;
					$prev_loc = $tokso[2];
					$start = 1;
				}
			}
			else{
				if($toksp[$pind2] == $tokso[$pindoff]){
					print OUTFILE "$gp_pat" . "_" ."$paternal\t";
					print OUTFILE "$tokso[2]\t";
					$currg = $gp_pat;
					$prev_loc = $tokso[2];
					$start = 1;
				}
			}

		}#end while loop for first informative SNP

		$switch_now = 0;
		while($linep = <INFILE_P>){
			chop $linep;
			@toksp = split('\t', $linep);
			$lineo = <INFILE_O>;
			chop $lineo;
			@tokso = split('\t', $lineo);

			if(($linep =~ /[emu]/) || ($lineo =~ /[emu]/)){				
				next;
			}
			if($toksp[$pind1] == $tokso[$pindoff]){				
				if($toksp[$pind2] == $tokso[$pindoff]){											
					next;
				}
				else{
					if($currg eq $gm_pat){						
						$prev_prev_loc = $prev_loc;
						$prev_loc = $tokso[2];
						$switch_now =0;
						next;
					}
					else{	
						if($switch_now == 0){
							$end_loc = $prev_loc;
							$start_loc = $tokso[2];
						}
						#print STDOUT "Entered switch: current kid: $tmp_ary[$i+1]\t current parent: $tmp_ary[$i]\t current gp: $gm_pat \t switch now =$switch_now\n";	
						$switch_now++;	
						$prev_loc = $tokso[2];
						$prev_prev_loc = $prev_loc;
						if($switch_now == $real_size){	
							#print STDOUT "Made switch: switch now =$switch_now\n";				
							print OUTFILE "$end_loc\n";
							print OUTFILE "$gm_pat" . "_" ."$paternal\t";
							print OUTFILE "$start_loc\t";
							$currg = $gm_pat;
							$switch_now = 0;
						}
					}
				}
			}
			else{
				if($toksp[$pind2] == $tokso[$pindoff]){					
					if($currg eq $gp_pat){
						$prev_prev_loc = $prev_loc;
						$prev_loc = $tokso[2];
						next;
						$switch_now = 0;
					}
					else{
						if($switch_now == 0){
							$end_loc = $prev_loc;
							$start_loc = $tokso[2];
						}
						#print STDOUT "Entered switch: current kid: $tmp_ary[$i+1]\t current parent: $tmp_ary[$i]\t current gp: $gp_pat \t switch now =$switch_now\n";
						$switch_now++;
						$prev_loc = $tokso[2];
						$prev_prev_loc = $prev_loc;
						if($switch_now == $real_size){
							#print STDOUT "Made switch: switch now =$switch_now\n";
							print OUTFILE "$end_loc\n";
							print OUTFILE "$gp_pat" . "_" ."$paternal\t";
							print OUTFILE "$start_loc\t";
							$currg = $gp_pat;
							$switch_now = 0;
						}
					}
				}				
			}

		}
	   print OUTFILE "$prev_loc\n";
	   close INFILE_P;
	   close INFILE_O;
	   close OUTFILE;
	   }
	   
	}
}




#write 3_gen_crossover file for maternal
foreach $offspr (keys %offspr_hash){
	if(exists $offspr_mat_hash{$offspr}){
		#write 3_gen_crossover file for maternal part
		@tmp_ary = @{$offspr_mat_hash{$offspr}};


	   $numtrios = @tmp_ary;

	   for($i = 0; $i < $numtrios; $i = $i+2){
		$trioline = $tmp_ary[$i+1];
		chop $line;
		@triotoks = split('_', $trioline);
		$maternal = $triotoks[0];
		$offspring = $triotoks[2];

		$fname_out = "$offspring" . "_crossovers_3gen.txt";
		if(-e $fname_out){
			open (OUTFILE, ">>$fname_out") || die;
		}
		else{
			open (OUTFILE, ">$fname_out") || die;
		}

		#get maternal's mat and pat (gp_mat, gm_mat)
		$line = $tmp_ary[$i];
		chop $line;
		@toks = split('_', $line);
		$gm_mat = $toks[0];
		$gp_mat = $toks[1];




		$fname_mat = "$tmp_ary[$i]" . "_geno.txt";
		$fname_offspr = "$tmp_ary[$i+1]" . "_geno.txt";
		
		open (INFILE_M, "$fname_mat") || die;
		open (INFILE_O, "$fname_offspr") || die;

		#traverse to first crossover unambiguous chr
		#read in header line from both files
		$line = <INFILE_M>; $line = <INFILE_O>;

		#set index for inherited from gmat or gpat for middle paternal
		$mind1 = 6;
		$mind2 = 8;
		$pindoff = 6;

		$currg = "";

		#boolean flag to first informative SNP
		$start = 0;
		$switch_now = 0;
		while($start == 0){
			$linem = <INFILE_M>;
			chop $linem;
			@toksm = split('\t', $linem);
			$lineo = <INFILE_O>;
			chop $lineo;
			@tokso = split('\t', $lineo);

			if(($linem =~ /[emu]/) || ($lineo =~ /[emu]/)){				
				next;
			}

			if($toksm[$mind1] == $tokso[$pindoff]){
				if($toksm[$mind2] == $tokso[$pindoff]){				
					next;
				}
				else{
					print OUTFILE "$gm_mat" . "_" ."$maternal\t";
					print OUTFILE "$tokso[2]\t";
					$currg = $gm_mat;
					$prev_loc = $tokso[2];
					$start = 1;
				}
			}
			else{
				if($toksm[$mind2] == $tokso[$pindoff]){
					print OUTFILE "$gp_mat" . "_" ."$maternal\t";
					print OUTFILE "$tokso[2]\t";
					$currg = $gp_mat;
					$prev_loc = $tokso[2];
					$start = 1;
				}
			}

		}#end while loop for first informative SNP

		while($linem = <INFILE_M>){
			chop $linem;
			@toksm = split('\t', $linem);
			$lineo = <INFILE_O>;
			chop $lineo;
			@tokso = split('\t', $lineo);

			if(($linem =~ /[emu]/) || ($lineo =~ /[emu]/)){				
				next;
			}
			if($toksm[$mind1] == $tokso[$pindoff]){				
				if($toksm[$mind2] == $tokso[$pindoff]){											
					next;
				}
				else{
					if($currg eq $gm_mat){						
						$prev_prev_loc = $prev_loc;
						$prev_loc = $tokso[2];
						$switch_now = 0;
						next;
					}
					else{	
						if($switch_now == 0){
							$end_loc = $prev_loc;
							$start_loc = $tokso[2];
						}
						#print STDOUT "Entered switch: current kid: $tmp_ary[$i+1]\t current parent: $tmp_ary[$i]\t current gp: $gm_mat \t switch now =$switch_now\n";
						$switch_now++;
						$prev_loc = $tokso[2];
						$prev_prev_loc = $prev_loc;
						if($switch_now == $real_size){
							#print STDOUT "Made switch: switch now =$switch_now\n";			
							print OUTFILE "$end_loc\n";
							print OUTFILE "$gm_mat" . "_" ."$maternal\t";
							print OUTFILE "$start_loc\t";
							$currg = $gm_mat;
							$switch_now = 0;
						}
						
					}
				}
			}
			else{
				if($toksm[$mind2] == $tokso[$pindoff]){					
					if($currg eq $gp_mat){
						$prev_prev_loc = $prev_loc;
						$prev_loc = $tokso[2];
						$switch_now = 0;
						next;
					}
					else{
						if($switch_now == 0){
							$end_loc = $prev_loc;
							$start_loc = $tokso[2];
						}
						#print STDOUT "Entered switch: current kid: $tmp_ary[$i+1]\t current parent: $tmp_ary[$i]\t current gp: $gp_pat \t switch now =$switch_now\n";
						$switch_now++;
						$prev_loc = $tokso[2];
						$prev_prev_loc = $prev_loc;
						if($switch_now == $real_size){
							#print STDOUT "Made switch: switch now =$switch_now\n";
							print OUTFILE "$end_loc\n";
							print OUTFILE "$gp_mat" . "_" ."$maternal\t";
							print OUTFILE "$start_loc\t";
							$currg = $gp_mat;
							$switch_now = 0;
						}
					}
				}				
			}

		}
	   print OUTFILE "$prev_loc\n";
	   close INFILE_M;
	   close INFILE_O;
	   close OUTFILE;

	   }
	   
	}

}