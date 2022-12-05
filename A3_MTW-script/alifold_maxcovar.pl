#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $maxcovar=0;

while(<>){
  next unless (/^\s+\d+/);
  my @data = split;
  next if ($data[@data-1] =~ /[\+\-]/);
#  print "-- $data[2]\n";
 if ($data[2]<=2){ # mismatches
  if (@data-5>$maxcovar){
    $maxcovar = @data-5;
  }
#  print join ("\t",@data),"\t",$#data,"\n";
  }
}
#  print "==> maxcovar $maxcovar\n";
exit($maxcovar);
