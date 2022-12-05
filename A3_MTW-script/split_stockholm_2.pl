#!/usr/bin/env perl
# Last changed Time-stamp: <2018-07-27 22:38:46 mtw>
# -*-CPerl-*-
#
# usage: split_stockholm.pl -a myln.stk
#

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Path::Class;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $show_version = 0;
my $VERSION="0.1";
my ($id,$alifile);
my $nr = 1;
my $accession=undef;
my %stks=();

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1)
  unless GetOptions("aln|a=s"    => \$alifile,
		    "accession"  => sub{$accession=1},
		    "version"    => sub{$show_version = 1},
		    "man"        => sub{pod2usage(-verbose => 2)},
		    "help|h"     => sub{pod2usage(1)}
		   );

if ($show_version == 1){
  print "split_stockholm $VERSION\n";
  exit(0);
}

unless (-f $alifile){
  warn "Could not find input file provided via --aln|-a option";
  pod2usage(-verbose => 0);
}

open my $input, "<", $alifile or die $!;
my $fn = "$nr.aln";

# parse multi-Stockholm file into %stks hash
while(<$input>){
  chomp;
  my $line = $_;
  unless ($line eq "//"){
    if ($line =~ /\#=GF\sAC\s+(.+)/){
      $accession = $1;
    }
    if ($line =~ /\#=GF\sID\s+(.+)/){
      $id = $1;
 #     print "+++++++++++ $id +++++++++++\n";
    }
    push @{$stks{$nr}{data}}, $line;
    if (defined $accession) {$id=$accession}
    $stks{$nr}{id}=$id;
    #print "$line\n";
  }
  else{
    $nr++;
 #   print "---------------------------------\n";
  }
}
close $input;

# process %stks hash and write individual Stockholm files
foreach my $i (keys %stks){
  my $name = $stks{$i}{id}.".stk";
  open my $alnfile, ">", $name or die $!;
  foreach my $l ( @{$stks{$i}->{data}} ){
    print $alnfile $l."\n";
  }
  print $alnfile "//\n";
  #print Dumper (\@{$stks{$i}->{data}});
  close $alnfile;
  #print "----\n";
}

__END__

=head1 NAME

split_stockholm.pl - Split multi-Stockholm sequence alignments

=head1 SYNOPSIS

split_stockholm.pl [--aln|-a I<FILE>]

=head1 DESCRIPTION

This tool splits a multi-Stockholm multiple sequence aligbnment (MSA)
file into individual Stockholm files that contain only a single MSA.

=head1 OPTIONS

=over

=item B<--aln|-a>

Multi-Stockholm MSA file

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt> and
E<lt>michael.wolfinger@univie.ac.atE<gt>

=cut
(END)
