#!/usr/local/bin/perl

### program finds founders, stores the triofile of these founders

#ARGV[0] is pedigree number
$ped = $ARGV[0];

#open list of all trios that have the offspring phased, store in @trios_ary
$filename_cleared = "Cleared_Ped_" . "$ped" . ".txt";
open (INFILE_CLEARED, $filename_cleared) || die;	

#open founder outfile
open (OUTFILE, ">GAW_nonfounder_matrix.txt") || die;

@trios_ary = <INFILE_CLEARED>;
$filelen = @trios_ary;

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


#find all children of each founder id's by traversing @trios_ary and storing trio file in @founder_kids

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


#print out column headers 
#TODO: ADD CHR 
print OUTFILE "Individ\tFounding_chr\tStart\tEnd\tPath\n";



foreach $founder (keys %second_gen_hash){
	@tmp_ary = @{$second_gen_hash{$founder}};
	$len = @tmp_ary;

	for($i = 0; $i < $len; $i++){
		@triotoks = split('_', $tmp_ary[$i]);
		$maternal = @triotoks[0];
		$paternal = @triotoks[1];
		$offspr = @triotoks[2];

		if($len >2){	
			$fnamecross = "$offspr" . "_crossover.txt";
			if (-e $fnamecross){
				open (INFILE_CROSS, "$fnamecross");
				while($line = <INFILE_CROSS>){
					chop $line;
					@toks = split("\t", $line);
					if ($toks[0] =~ /$founder/){
						print OUTFILE "$offspr\t";
						print OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
						$path = "$offspr" . "_" . "$toks[0]";
						print OUTFILE "$path\n";
					}
				}
				close INFILE_CROSS;
			}


						
		}
		elsif($len == 2){
			$fnamecross = "$offspr" . "_crossover.txt";
			if (-e $fnamecross){
				open (INFILE_CROSS, "$fnamecross");
				while($line = <INFILE_CROSS>){
					chop $line;
					@toks = split("\t", $line);
					if ($toks[0] =~ /$founder/){
						print OUTFILE "$offspr\t";
						print OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
						$path = "$offspr" . "_" . "$toks[0]";
						print OUTFILE "$path\n";
					}
				}
				close INFILE_CROSS;
			}
		}
		else{
			#this founder only has either one or two children since no crossover file exists
			$fname_geno = "$tmp_ary[0]" . "_geno.txt";
			open (INFILE_GENO, $fname_geno);
			$headerline = <INFILE_GENO>;
			@lines = <INFILE_GENO>;
			$num_lines = @lines;						

			print OUTFILE "$offspr\t";
			print OUTFILE "$founder" . "_" . "1\t";
			@toks = split("\t", $lines[0]);
			$start = $toks[2];
			print OUTFILE "$start\t";
			@toks = split("\t", $lines[$num_lines-1]);
			$end = $toks[2];
			print OUTFILE "$end\t";
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
	my $fname = "$this_offspr" . "_crossovers_3gen.txt";
	open (INFILETHREE, $fname) || die;
	
	open (OUTFILE, ">>GAW_nonfounder_matrix.txt") || die;	

	my @lines3gen = <INFILETHREE>;
	my $line3;
	my $i;
	my $len = @lines3gen;

	for($i = 0; $i < $len; $i++){
		$line3 = $lines3gen[$i];
		chop $line3;
		my @toks = split('\t', $line3);
		my $pair3 = $toks[0];
		@toks2 = split('_', $toks[0]);

		#rearranging order of 3gen file ids to match nonfounder matrix
		my $id1 = $toks2[1];
		my $id2 = $toks2[0];

		open (INFILEMATRIX, "GAW_nonfounder_matrix.txt") || die;
		

		my $linemat = <INFILEMATRIX>;

		while($linemat = <INFILEMATRIX>){
			chop $linemat;
			my @toksmat = split('\t', $linemat);
			@toksmat2 = split('_', $toksmat[4]);
			my $id1mat = $toksmat2[0];
			my $id2mat = $toksmat2[1];
			my @tmp_ary;

			my $b3 = $toks[1];
			my $c3 = $toks[2];
			my $d3 = $toksmat[2];
			my $e3 = $toksmat[3];			

			my $newpath = "$this_offspr" . "_" . "$toksmat[4]";

			#print "$this_offspr\t$b3\t$c3\t$d3\t$e3\t$newpath\n";

			

			if(($id1 eq $id1mat) && ($id2 eq $id2mat)){
				if($b3 > $e3){
					next;
				}
				elsif($c3 < $d3){
					next;
				}
				elsif($b3 < $d3){
					if($c3 == $d3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$d3\t$d3\t";
						print OUTFILE "$newpath\n";
					}
					elsif($c3 < $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$d3\t$c3\t";
						print OUTFILE "$newpath\n";

					}
					elsif($c3 == $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$d3\t$c3\t";
						print OUTFILE "$newpath\n";

					}
					elsif($c3 > $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$d3\t$e3\t";
						print OUTFILE "$newpath\n";

					}
					else{
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "Error\tError\t";
						print OUTFILE "$newpath\n";
					}
				}
				elsif($b3 == $d3){
					if($c3 < $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$b3\t$c3\t";
						print OUTFILE "$newpath\n";

					}
					elsif($c3 == $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$b3\t$c3\t";
						print OUTFILE "$newpath\n";

					}
					elsif($c3 > $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$b3\t$e3\t";
						print OUTFILE "$newpath\n";

					}
					else{
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "Error\tError\t";
						print OUTFILE "$newpath\n";

					}
				}
				elsif($b3 > $d3){
					if($c3 < $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$b3\t$c3\t";
						print OUTFILE "$newpath\n";

					}
					elsif($c3 == $e3){
						print OUTFILE "$this_offspr\t$toksmat[1]\t";
						print OUTFILE "$b3\t$c3\t";
						print OUTFILE "$newpath\n";

					}
					elsif($c3 > $e3){
						if($b3 == $e3){
							print OUTFILE "$this_offspr\t$toksmat[1]\t";
							print OUTFILE "$b3\t$e3\t";
							print OUTFILE "$newpath\n";

						}
						else{
							print OUTFILE "$this_offspr\t$toksmat[1]\t";
							print OUTFILE "$b3\t$e3\t";
							print OUTFILE "$newpath\n";

						}
					}
				}
				else{
					print OUTFILE "$this_offspr\t$toksmat[1]\t";
					print OUTFILE "ERROR\t$newpath\n";
				}
							
			}			

		}

		close INFILEMATRIX;
		
		
	}

	close INFILETHREE;	
	close OUTFILE;
}




















