#!/usr/local/bin/perl


open (OUTFILE, ">out_012.txt") || die;

$pednumber = $ARGV[0];
$chrnumber = $ARGV[1];

$INLOCSNAME = "chr"."$chrnumber"."_ped"."$pednumber".".012.pos";
$INGENONAME = "chr"."$chrnumber"."_ped"."$pednumber".".012_transposed";
$INIDNAME = "chr"."$chrnumber"."_ped"."$pednumber".".012.indv";

open (INLOCS, $INLOCSNAME) || die;
open (INGENO, $INGENONAME) || die;
open (INID, $INIDNAME) || die;


print OUTFILE "Name\tChr\tPosition\t";

while($line = <INID>){
	chop $line;
	print OUTFILE "$line\t";
}
print OUTFILE "\n";

close INID;

#read in header line
$line = <INGENO>;

while($lineloc = <INLOCS>){
	chop $lineloc;
	@toklocs = split('\t', $lineloc);
	$marker = "$toklocs[0]" . "_" . "$toklocs[1]";

	print OUTFILE "$marker\t";

	print OUTFILE "$toklocs[0]\t$toklocs[1]\t";

	$linegeno = <INGENO>;
	$linegeno =~ s/-1/-/g;

	print OUTFILE $linegeno;
}

close INLOCS;
close INGENO;
close OUTFILE;













