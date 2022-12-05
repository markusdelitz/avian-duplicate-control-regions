#!/usr/bin/env perl
# Last changed Time-stamp: <2020-06-21 23:45:32 mtw>
# strip MSA, i.e. remove redundant sequences

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Carp;
use Data::Dumper;
use Pod::Usage;
use Bio::AlignIO;
use diagnostics;

my %format = (
              'C' => 'clustalw',
              'S' => 'stockholm',
              'M' => 'maf',
              'F' => 'fasta',
             );

my $infile_aln = undef;
my $alnformat = undef;
my ($in,$out,$aln);
my @keep = ();
my %seen = ();

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1)
 unless GetOptions("a|aln=s"    => \$infile_aln,
                   "f|format=s" => \$alnformat,
                   "man"        => sub{pod2usage(-verbose => 2)},
                   "help|h"     => sub{pod2usage(1)}
                   );

unless (-f $infile_aln){
  warn "Cannot find input alignment privided via -a|--aln option ...";
  pod2usage(-verbose => 0);
}
unless (defined $alnformat){
  warn "Format of alignment not given; Please use -f|--format [C|S|M|F] option";
  pod2usage(-verbose => 0);
}
if (defined $format{$alnformat}){
  $alnformat = $format{$alnformat};
}
else {
  croak "ERROR: input alignment format $alnformat not supported.
Accepted formats are 'C' (Clustal), 'S' (Stockholm), 'M' (Mafft) or 'F' (Fasta).";
}

if (defined $infile_aln){
  $in  = Bio::AlignIO->new(-file   => "$infile_aln",
                           -format => $alnformat,
                           -displayname_flat => 1);
  $out = Bio::AlignIO->new(-fh   => \*STDOUT ,
                           -format => $alnformat);
}
else { croak "Could not process input alignment" }

$aln = $in->next_aln();
#print Dumper($aln);
$aln->set_displayname_flat();

foreach my $i (1..$aln->num_sequences){
  my $alnT = $aln->get_seq_by_pos($i);
  my $seq = $alnT->seq();
  if (exists($seen{$seq})){next}
  else {
  push @keep, $i;
#  print "keep $seq\n";
  }
  $seen{$seq}=1;
}

my $subset = $aln->select_noncont(@keep);
$out->write_aln($subset);

__END__

=head1 NAME

strip_aln.pl - Remove redundant sequences from multiple sequence alignments

=head1 SYNOPSIS

strip_aln.pl [-a|--aln I<FILE>]  [options]

=head1 DESCRIPTION

This script removes redundant lines from a multiple sequence alignment (MSA).
The primary purpose of this is to obtain stripped MSAs for follow-up
consensus structure prediction or covariance model computation.

=head1 OPTIONS

=over

=item B<-a|--aln>

Input file in Clustal format. Either this option or B<-f|--fa> mus be
given.


=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
