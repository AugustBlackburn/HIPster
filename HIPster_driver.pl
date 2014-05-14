#!/usr/bin/perl

use HIPster
$data=$ARGV[0];
$pedigree=$ARGV[1];
$markers=$ARGV[2];

&Step1($data,$pedigree);
&Step2($markers);
&Step3($markers);
&Step4($markers);
&Step5();
&Step6();
&Step7($data);
&Step8($data);
&Step9();
&Step10();
&Step11();
&Step12($data);
print STDOUT "Finished. Thanks for using HIPster!\n";

