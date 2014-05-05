#!/usr/local/bin/perl

### program phases both parents into two separate files named "[parentid]_phased.txt" and
### writes a file for each offspring of the ranges between each crossover event named "[offsprid]_crossover.txt"


#ARGV[0] is pedigree number
$ped = $ARGV[0];

#ARGV[1] is marker requirement for switches
$requirement = $ARGV[1];

if($requirement == 0){
	print STDOUT "A minimum of 1 marker is required to identify a crossover.";
	die;
}
if($requirement > 10){
	print STDOUT "You are using >10 markers to identify a crossover. If you observe strange results try lowering the number of markers.";
}

#set mother and father indices from "_geno"
$mind1 = 6;
$mind2 = 7;
$pind1 = 8;
$pind2 = 9;

#open list of phased trios and store data in file
$filename_cleared = "Cleared_Ped_" . "$ped" . ".txt";
open (INFILE_CLEARED, $filename_cleared) || die;
$tiolinecount = 0;
while($line = <INFILE_CLEARED>){
	chop $line;
	@triotoks = split('_', $line);
	@maternal[$triolinecount] = @triotoks[0];
	@paternal[$triolinecount] = @triotoks[1];
	@offspr[$triolinecount] = @triotoks[2];
	$triolinecount++;
}

#find all unique mothers
my @uniquemothers;
my %seenmothers;

foreach my $value (@maternal){
	if (! $seenmothers{$value}++ ){
		push @uniquemothers, $value;
	}
}

#find all unique fathers
my @uniquefathers;
my %seenfathers;

foreach my $value2 (@paternal){
	if (! $seenfathers{$value2}++ ){
		push @uniquefathers, $value2;
	}
}

#for each mother, apply minimum recombination model
foreach $currentmom (@uniquemothers){
	#create a list of the current mom's kids
	#print STDOUT "Current mom is: $currentmom\nHer kids are:\n";
	$numberkids = 0;
	for($i=0;$i<$triolinecount;$i++){
		if(@maternal[$i] eq $currentmom){
			@curroffspring[$numberkids]= @offspr[$i];
			@currdad[$numberkids]= @paternal[$i];
			$numberkids++;
			#print STDOUT "@offspr[$i]\n";
		}
	}
	#print STDOUT "She has $numberkids kids.\n\n";
	
	#if the mom doesn't have 2 kids move on
	if($numberkids != 2){
		next;
	}

	#open outfile for phased mom file
	$mfname = "$currentmom" . "_phased.txt";
	open (M_OUTFILE, ">$mfname") || die;
	
	@offspringinfilearray = ();
	#open infile and outfiles for kids
	for ($i = 0;$i<$numberkids;$i++){
		$kidinfilename = "$currentmom"."_"."@currdad[$i]"."_"."@curroffspring[$i]"."_geno.txt";
		$kidoutfilename = "@curroffspring[$i]"."_crossover.txt";
		#print STDOUT "$kidinfilename\n";
		open ($offspringinfilearray[$i], "<$kidinfilename") || die;
		open ($offspringoutfilearray[$i], ">$kidoutfilename") || die;	
	}
	print M_OUTFILE "SNP Name\t\Chr\tPosition\t";
	print M_OUTFILE "$currentmom" . "_1\t" . "$currentmom" . "_2\n";

	#boolean flag because homozygous maternal SNPs provide no information, need to find first heterozygous SNP
	$start = 0;

	#traverse to first heterozygous maternal SNP
	#read in header line 
	for($i = 0; $i < $numberkids; $i++){	
		$lines[$i] = readline($offspringinfilearray[$i]);
	}

	while ($start == 0){
		$minct = 0;
		for($i = 0; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chop $lines[$i];		
		}
		@toks = split('\t', $lines[0]);

		if($toks[3] != 1){
			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
			print M_OUTFILE "-\t-\n";
			next;
		}
		else{
			for($i = 0; $i < $numberkids; $i++){
				@toks = split('\t', $lines[$i]);
				if ($toks[$mind1] !~ /[emu]/){			
					$minct++;
				}
			}
			if($minct != $numberkids){
				print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
				print M_OUTFILE "-\t-\n";
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
		@toks = split('\t', $lines[$i]);
		if($toks[$mind1] == 1){
			#set to chr1 (index $mind1 in _geno file), else to chr2 (index $mind2)
			$chr_ary[$i] = "$mind1";
		}
		else{
			$chr_ary[$i] = "$mind2";
		}
	}
	print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t1\t0\n";
	#write out the parent chr_1 and beginning location of the first segment for each offspring file
	for($i = 0; $i < $numberkids; $i++){
		$prev_hetero_loc[$i] = $toks[2];
		$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];
		if($chr_ary[$i] == $mind1){
			print {$offspringoutfilearray[$i]} "$currentmom" . "_1\t" ."$toks[2]\t";				
		}
		else{
			print {$offspringoutfilearray[$i]} "$currentmom" . "_2\t" ."$toks[2]\t";				
		}
	}			
	#begin phasing after initializing parent and offspr's @chr_ary for each offspr's current chr1 or chr2
	while(!(eof $offspringinfilearray[0])){
		$minct = 0;
		
		#if parent is homozygous, we already know the output
		$lines[0] = readline($offspringinfilearray[0]);
		chop $lines[0];
		@toks = split('\t', $lines[0]);	
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

			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
			print M_OUTFILE "$toks[3]\t$toks[3]\n";
			#advance to next line for rest of offspr
			for($i = 1; $i < $numberkids; $i++){	
				$lines[$i] = readline($offspringinfilearray[$i]);
			}
			next;			
		}
		#parent is heterozygous
		#already read in line from first offspr so start at $i = 1, take care of first offspr first
		if ($toks[$mind1] !~ /[emu]/){							
			$val_ary[0] = $toks[$chr_ary[0]];
			$minct++;				
		}
		else{
			$val_ary[0] = "";
		}
		for($i = 1; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray[$i]);
			chop $lines[$i];
			@toks = split('\t', $lines[$i]);
			if ($toks[$mind1] !~ /[emu]/){							
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
				#print STDOUT "ct is $ct\tct2 is $ct2\n";
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
						#print {$offspringoutfilearray[$i]} "Entered common allele is 0 counts = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray[$i]} "Entered common allele is 0 switch count 1\n";
						}
						if($requirement==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							#print {$offspringoutfilearray[$i]} "Entered common allele is 0 switch\n";
							if($chr_ary[$i] == $mind1){
								$chr_ary[$i] = $mind2;
							}
							elsif($chr_ary[$i] == $mind2){
								$chr_ary[$i] = $mind1;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == $mind1){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1\t";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i]\t";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2\t";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i]\t";
							}	
						}
					}
				}
				elsif($ct2>$ct){
					#print STDOUT "Entered the commom allele is 1\n.";
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
						#print {$offspringoutfilearray[$i]} "Entered common allele is 1 counts = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray[$i]} "Entered common allele is 1 switch count 1\n";
						}
						if($requirement==$currcount1){
							@counts[$i] = ();
							#change crossover offspr strand	
							#print {$offspringoutfilearray[$i]} "Entered common allele is 1 switch\n";
							if($chr_ary[$i] == $mind1){
								$chr_ary[$i] = $mind2;
							}
							elsif($chr_ary[$i] == $mind2){
								$chr_ary[$i] = $mind1;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == $mind1){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1\t";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i]\t";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2\t";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i]\t";
							}		
						}
					}
				}
				elsif($ct2==$ct){
					#print STDOUT "Entered parsimony.\n";
					#throw out SNP, its ambiguous, what if this happens?
					#synthetic switch?
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
						if($requirement==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							#print {$offspringoutfilearray[$i]} "Entered parsimony switch\n";
							if($chr_ary[$i] == $mind1){
								$chr_ary[$i] = $mind2;
							}
							elsif($chr_ary[$i] == $mind2){
								$chr_ary[$i] = $mind1;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == $mind1){
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_1\t";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i]\t";
							}
							else{
								print {$offspringoutfilearray[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray[$i]} "$currentmom" . "_2\t";
								print {$offspringoutfilearray[$i]} "$newswitchbeginning[$i]\t";
							}	
						}
					}
				}
			}
			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t$m_chr1\t$m_chr2\n";	
		}
		#there weren't at least three offspr with genotype data
		else{
			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t-\t-\n";	
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
foreach $currentdad (@uniquefathers){
	#create a list of the current mom's kids
	#print STDOUT "Current dad is: $currentdad\nHis kids are:\n";
	$numberkids = 0;
	for($i=0;$i<$triolinecount;$i++){
		if(@paternal[$i] eq $currentdad){
			@curroffspring2[$numberkids]= @offspr[$i];
			@currmom[$numberkids]= @maternal[$i];
			$numberkids++;
			#print STDOUT "@offspr[$i]\n";
		}
	}
	#print STDOUT "He has $numberkids kids.\n\n";
	
	#if the dad doesn't have at least 3 kids move on
	if($numberkids != 2){
		next;
	}

	#open outfile for phased mom file
	$mfname = "$currentdad" . "_phased.txt";
	open (M_OUTFILE, ">$mfname") || die;
	
	@offspringinfilearray2 = ();
	#open infile and outfiles for kids
	for ($i = 0;$i<$numberkids;$i++){
		$kidinfilename = "@currmom[$i]"."_"."$currentdad"."_"."@curroffspring2[$i]"."_geno.txt";
		$kidoutfilename = "@curroffspring2[$i]"."_crossover.txt";
		#print STDOUT "$kidinfilename\n";
		open ($offspringinfilearray2[$i], "<$kidinfilename") || die;
		open ($offspringoutfilearray2[$i],">>", "$kidoutfilename") || die;	
	}
	#for ($i = 0;$i<$numberkids;$i++){
	#	print {$offspringoutfilearray2[$i]} "Can print to file.\n";
	#}
	print M_OUTFILE "SNP Name\t\Chr\tPosition\t";
	print M_OUTFILE "$currentdad" . "_1\t" . "$currentdad" . "_2\n";

	#boolean flag because homozygous maternal SNPs provide no information, need to find first heterozygous SNP
	$start = 0;

	#traverse to first heterozygous maternal SNP
	#read in header line 
	for($i = 0; $i < $numberkids; $i++){	
		$lines[$i] = readline($offspringinfilearray2[$i]);
	}

	while ($start == 0){
		$minct = 0;
		for($i = 0; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray2[$i]);
			chop $lines[$i];
			#print STDOUT "$lines[$i]\n";		
		}
		@toks = split('\t', $lines[0]);

		if($toks[4] != 1){
			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
			print M_OUTFILE "-\t-\n";
			next;
		}
		else{
			for($i = 0; $i < $numberkids; $i++){
				@toks = split('\t', $lines[$i]);
				if ($toks[$pind1] !~ /[emu]/){			
					$minct++;
				}
			}
			if($minct != $numberkids){
				print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
				print M_OUTFILE "-\t-\n";
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
	#print STDOUT "Phasing $currentdad.\n";
	for($i = 0; $i < $numberkids; $i++){
		@toks = split('\t', $lines[$i]);
		if($toks[$pind1] == 1){
			#set to chr1 (index $pind1 in _geno file), else to chr2 (index $pind2)
			$chr_ary[$i] = "$pind1";
		}
		else{
			$chr_ary[$i] = "$pind2";
		}
	}
	print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t1\t0\n";
	#write out the parent chr_1 and beginning location of the first segment for each offspring file
	for($i = 0; $i < $numberkids; $i++){
		$prev_hetero_loc[$i] = $toks[2];
		$prev_prev_hetero_loc[$i] = $prev_hetero_loc[$i];
		if($chr_ary[$i] == $pind1){
			print {$offspringoutfilearray2[$i]} "$currentdad" . "_1\t" ."$toks[2]\t";				
		}
		else{
			print {$offspringoutfilearray2[$i]} "$currentdad" . "_2\t" ."$toks[2]\t";				
		}
	}			
	#begin phasing after initializing parent and offspr's @chr_ary for each offspr's current chr1 or chr2
	while(!(eof $offspringinfilearray2[0])){
		$minct = 0;
		
		#if parent is homozygous, we already know the output
		$lines[0] = readline($offspringinfilearray2[0]);
		chop $lines[0];
		@toks = split('\t', $lines[0]);	
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

			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t";
			print M_OUTFILE "$toks[4]\t$toks[4]\n";
			#advance to next line for rest of offspr
			for($i = 1; $i < $numberkids; $i++){	
				$lines[$i] = readline($offspringinfilearray2[$i]);
			}
			next;			
		}
		#parent is heterozygous
		#already read in line from first offspr so start at $i = 1, take care of first offspr first
		if ($toks[$pind1] !~ /[emu]/){							
			$val_ary[0] = $toks[$chr_ary[0]];
			$minct++;				
		}
		else{
			$val_ary[0] = "";
		}
		for($i = 1; $i < $numberkids; $i++){	
			$lines[$i] = readline($offspringinfilearray2[$i]);
			chop $lines[$i];
			@toks = split('\t', $lines[$i]);
			if ($toks[$pind1] !~ /[emu]/){							
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
				#print STDOUT "ct is $ct\tct2 is $ct2\n";
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
						#print {$offspringoutfilearray2[$i]} "Entered common allele is 0 counts = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray2[$i]} "Entered common allele is 0 switch count 1\n";
						}
						if($requirement==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							#print {$offspringoutfilearray2[$i]} "Entered common allele is 0 switch\n";
							if($chr_ary[$i] == $pind1){
								$chr_ary[$i] = $pind2;
							}
							elsif($chr_ary[$i] == $pind2){
								$chr_ary[$i] = $pind1;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == $pind1){
								print {$offspringoutfilearray2[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray2[$i]} "$currentdad" . "_1\t";
								print {$offspringoutfilearray2[$i]} "$newswitchbeginning[$i]\t";
							}
							else{
								print {$offspringoutfilearray2[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray2[$i]} "$currentdad" . "_2\t";
								print {$offspringoutfilearray2[$i]} "$newswitchbeginning[$i]\t";
							}	
						}
					}
				}
				elsif($ct2>$ct){
					#print STDOUT "Entered the commom allele is 1\n.";
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
						#print {$offspringoutfilearray2[$i]} "Entered common allele is 1 counts = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray2[$i]} "Entered common allele is 1 switch count 1\n";
						}
						if($requirement==$currcount1){
							@counts[$i] = ();
							#change crossover offspr strand	
							#print {$offspringoutfilearray2[$i]} "Entered common allele is 1 switch\n";
							if($chr_ary[$i] == $pind1){
								$chr_ary[$i] = $pind2;
							}
							elsif($chr_ary[$i] == $pind2){
								$chr_ary[$i] = $pind1;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == $pind1){
								print {$offspringoutfilearray2[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray2[$i]} "$currentdad" . "_1\t";
								print {$offspringoutfilearray2[$i]} "$newswitchbeginning[$i]\t";
							}
							else{
								print {$offspringoutfilearray2[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray2[$i]} "$currentdad" . "_2\t";
								print {$offspringoutfilearray2[$i]} "$newswitchbeginning[$i]\t";
							}		
						}
					}
				}
				elsif($ct2==$ct){
					#print STDOUT "Entered parsimony.\n";
					#throw out SNP, its ambiguous, what if this happens?
					#synthetic switch?
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
						#print {$offspringoutfilearray2[$i]} "Entered common allele parsimonious = @counts[$i]\n";
						if($currcount1 == $currcountone){
							$oldswitchend[$i] = $prev_prev_hetero_loc[$i];
							$newswitchbeginning[$i] = $prev_hetero_loc[$i];
							#print {$offspringoutfilearray2[$i]} "Entered parsimonious switch count 1\n";
						}
						if($requirement==$currcount1){
							#change crossover offspr strand
							@counts[$i] = ();
							#print {$offspringoutfilearray2[$i]} "Entered parsimony switch\n";
							if($chr_ary[$i] == $pind1){
								$chr_ary[$i] = $pind2;
							}
							elsif($chr_ary[$i] == $pind2){
								$chr_ary[$i] = $pind1;
							}
							#output if on chr1 or chr2 of parent
							if($chr_ary[$i] == $pind1){
								print {$offspringoutfilearray2[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray2[$i]} "$currentdad" . "_1\t";
								print {$offspringoutfilearray2[$i]} "$newswitchbeginning[$i]\t";
							}
							else{
								print {$offspringoutfilearray2[$i]} "$oldswitchend[$i]\n";
								print {$offspringoutfilearray2[$i]} "$currentdad" . "_2\t";
								print {$offspringoutfilearray2[$i]} "$newswitchbeginning[$i]\t";
							}	
						}
					}
				}
			}
			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t$m_chr1\t$m_chr2\n";	
		}
		#there weren't at least three offspr with genotype data
		else{
			print M_OUTFILE "$toks[0]\t$toks[1]\t$toks[2]\t-\t-\n";	
		}
	}
	#print end of last segment for each offspring (last heterozygous SNP and close the files
	for($i = 0; $i < $numberkids; $i++){
		print {$offspringoutfilearray2[$i]} "$prev_hetero_loc[$i]\n";
		close $offspringoutfilearray2[$i];
	}
	for($i = 0; $i < $numberkids; $i++){
		close $offspringinfilearray2[$i];
	}
	close M_OUTFILE;
}













