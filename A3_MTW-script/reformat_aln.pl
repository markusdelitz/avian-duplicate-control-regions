#!/usr/bin/env perl
# Last changed Time-stamp: <2018-12-30 17:14:36 mtw>
# -*-CPerl-*-
#
# reformat_aln.pl: inter-convert alignment formats

use strict;
use warnings;
use Bio::AlignIO;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;

my $infile_aln = undef;
my $infmt = undef;
my $outfmt = undef;
my $display_flat = 0;
my $display_count = 0;

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless
    GetOptions("a|aln=s"        => \$infile_aln,
	       "i|infmt=s"      => \$infmt,
	       "o|outfmt=s"     => \$outfmt,
	       "f|displayflat"  => sub{$display_flat=1},
	       "c|displaycount" => sub{$display_count=1},
	       "man"        => sub{pod2usage(-verbose => 2)},
	       "help|h"     => sub{pod2usage(1)}
    );

unless (-f $infile_aln){
  warn "Could not find input alignment file provided via -a|--aln option";
  pod2usage(-verbose => 0);
}

if ($display_flat==1 && $display_count==1){$display_flat=0}

my $in = Bio::AlignIO->new(-file => "$infile_aln" ,
			   -format => $infmt);

my $out = Bio::AlignIO->new(-fh   => \*STDOUT ,
			    -format => $outfmt);

while ( my $aln = $in->next_aln ) {
  if ( $display_flat ) {
    $aln->set_displayname_flat(1);
  }
  if ( $display_count ) {
    $aln->set_displayname_count(1);
  }
  $out->write_aln($aln);
}

=head1 NAME

reformat_aln.pl - convert alignment file formats

=head1 SYNOPSIS

reformat_aln.pl [-a|--aln I<FILE>] [-i|--infmt I<STRING>] 
[-o|--outfmt I<STRING>] [options]

=head1 DESCRIPTION

This tool is a simple L<Bio::AlignIO> based converter for various
sequence alignment file formats. Among others, these file formats are available:

=over

=item fasta

=item mase

=item selex

=item clustalw

=item msf

=item phylip

=item po

=item stockholm

=item XMFA

=item metafasta

=back

See the L<AlignIO and SimpleAlign
HOWTO|https://bioperl.org/howtos/AlignIO_and_SimpleAlign_HOWTO.html>
for a full list of L<Bio::AlignIO> supported file formats.

=head1 OPTIONS

=over

=item B<--aln|-a>

Inpt alignment file

=item B<--infmt|-i>

Input file format

=item B<--outfmt|-o>

Output file format

=item  B<--displayflat|-f>

Makes all the sequences be displayed as just their name, not
name/start-end

=item  B<--displaycount|-c>

Sets the names to be name_\# where \# is the number of times this name
\# has been used.

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt> and
E<lt>michael.wolfinger@univie.ac.atE<gt>

=cut

