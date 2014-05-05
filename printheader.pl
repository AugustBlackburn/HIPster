#!/usr/bin/perl -w

$datafile = $ARGV[0];

open DATAFILE, $datafile or die;
open OUTFILE, ">headerline.txt" or die;

$dataheader = <DATAFILE>;
print OUTFILE "$dataheader\n";

close OUTFILE;