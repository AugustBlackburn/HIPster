#!/usr/local/bin/perl

### program finds diploid markers from founder & non-founder matrices

open (INFILE_F, "imputed_founder_chrs.txt");
open (INFILE_N, "GAW_nonfounder_matrix.txt");

open (OUTFILE, ">GAW_diploid_matrix.txt");

#read in header line of founder matrix
$line = <INFILE_F>;
$len = length($line);

#get hash of indices for each founder id with this header line
%founder_ind = ();
@toks = split('\t', $line);

for($i = 3; $i < $len; $i++){
	$founder_ind{$toks[$i]} = $i;
}

#make hash of non-founders
%non_founders = ();

#get rid of non_founder header line
$line = <INFILE_N>;

while($line = <INFILE_N>){
	@toks = split('\t', $line);
	$indiv = $toks[0];
	$non_founders{$indiv} = 1;
}

#print OUTFILE header line
print OUTFILE "SNP\tChr\tBp\t";

foreach $non_f (keys %non_founders){
	print OUTFILE "$non_f\t";
#print "first" . "$non_f\t";
}
#print "\n";

print OUTFILE "\n";


close INFILE_N;
open (INFILE_N, "GAW_nonfounder_matrix.txt");
#get rid of non_founder header line
$line = <INFILE_N>;

@non_flines = <INFILE_N>;
close INFILE_N;

#output SNP, Chr and bp from founder matrix
while($linef = <INFILE_F>){
	@toksf = split('\t', $linef);
	print OUTFILE "$toksf[0]\t$toksf[1]\t$toksf[2]\t";
	$bp = $toksf[2];

	#index for output of each non_founder


	foreach $non_f (keys %non_founders){
#print "second" . "$non_f\t";

		open (INFILE_N, "GAW_nonfounder_matrix.txt");
		#read in header file
		$line = <INFILE_N>;

		$sum = 0;

		$founder1 = "";
		$founder2 = "";
		$has_f = 0;
		while($line = <INFILE_N>){
			
			@toks = split('\t', $line);
			if ($toks[0] eq $non_f){
				if(($bp >= $toks[2]) && ($bp <= $toks[3])){
					if($has_f == 0){
						$founder1 = $toks[1];
						$has_f = 1;
					}
					else{
						$founder2 = $toks[1];
						$has_f = 2;
						#print "$line\n";
					}
				}
			}
		}
		#There weren't 2 founders, no diploid entry
		if ($has_f < 2){
			print OUTFILE "-\t";
		}
		else{
			$ind1 = $founder_ind{$founder1};
			$ind2 = $founder_ind{$founder2};
	
			#print "$non_f\t$bp\t$founder1\t$founder2\n";
			#print "$founder_ind{$founder1}\t$founder_ind{$founder2}\n";




			if (($toksf[$ind1]=~ /[0-9]/) && ($toksf[$ind2]=~ /[0-9]/)) {
			

				$sum = $toksf[$ind1] + $toksf[$ind2];
				#print "$toksf[$ind1],$toksf[$ind2],$sum\n";
				print OUTFILE "$sum\t";
			}
			else{
				print OUTFILE "-\t";
			}
		}
#print "\n";

		close INFILE_N;
	}

	print OUTFILE "\n";
}
			
			
		
			
	















