#!/usr/bin/perl
 
use strict;
use warnings;
 
my $file = "GAW_nonfounder_matrix.txt";
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