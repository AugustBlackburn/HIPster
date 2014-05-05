#!/usr/bin/perl

# Transpose line <-> columns
# Didier Gonze
# 12/12/2003

$infile = "";
$outfile = "";

###  READ ARGUMENTS  ##################################################

foreach my $a (0..$#ARGV) {

   ### input directory 
   if ($ARGV[$a] eq "-i") {
      $infile = $ARGV[$a+1];
   }

   ### output file 
   elsif ($ARGV[$a] eq "-o") {
      $outfile = $ARGV[$a+1];
   }

   ### help
   elsif ($ARGV[$a] eq "-h") {
      die "Syntax: transpose.pl -i inputfilename -o outputfilename (default output=inputfilename.tr)\n";
   }

}

if ($infile eq ""){die "STOP! You have to specify the input file name";}
if ($outfile eq ""){$outfile="$infile".".tr"}


###  READ DATA  ##################################################

open inf, $infile or die "STOP! File $infile not found";
open ouf, ">$outfile";

$i=0;

foreach $line (<inf>){
  chomp $line;
  @line=split /\t/,$line;
  for ($j=0;$j<=$#line;$j++){
    $R[$i][$j]=$line[$j];
  }
  $i++;
}

$nlin=$i-1;
$ncol=$#line;

###  WRITE DATA  ##################################################

for ($i=0;$i<=$ncol;$i++){
  for ($j=0;$j<=$nlin;$j++){
    print ouf "$R[$j][$i]";
    if ($j<$nlin){
      print ouf "\t";
    }
  }
  print ouf "\n";
}

print "Output file $outfile created\n";

close inf;
close ouf;