#!/usr/bin/env perl
#
# remove gap-only sequences from  Clustal files

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Carp;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Bio::AlignIO;

my $infile_aln = undef;
my ($in,$out);

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("a|aln=s"    => \$infile_aln,
                                           "man"        => sub{pod2usage(-verbose => 2)},
                                           "help|h"     => sub{pod2usage(1)}
					  );

unless (defined ($infile_aln)){
  warn "Please provide an infput alignment in ClustalW format to  -a|--aln option ...";
  pod2usage(-verbose => 0);
}

unless (-f $infile_aln){
  warn "Cannot find input alignment privided via-a|--aln option ...";
  pod2usage(-verbose => 0);
}

if (defined $infile_aln){
  $in  = Bio::AlignIO->new(-file   => "$infile_aln",
			   -format => 'ClustalW',
			   -displayname_flat => 1);
  $out = Bio::AlignIO->new(-fh   => \*STDOUT ,
			   -format => 'ClustalW');
}
else { croak "Could not process input alignment" }

my $aln = $in->next_aln();
$aln->set_displayname_flat();
my $l   =  $aln->length;
my $dim =  $aln->num_sequences;
#print "+++ $dim sequences with length $l  in alignment +++\n";
my $gapstring = ("-"x$l);
#print "gapstring\n$gapstring\n\n\n";
my @keep = ();

foreach my $i (1..$dim){
  my $alnT = $aln->get_seq_by_pos($i);
  my $seq = $alnT->seq();
  # check if we have a gap-only substring of length $l
  index($seq,$gapstring) == 0 ? next : 1;
  push @keep, $i;
#  print "keep $seq\n";
}

my $subset = $aln->select_noncont(@keep);
$out->write_aln($subset);

__END__

=head1 NAME

remove-gaponly.pl - Remove gap-only sequences from ClustalW files

=head1 SYNOPSIS

remove-gaponly.pl [-a|--aln I<FILE>]  [options]

=head1 DESCRIPTION

This script removes gap-only lines from Clustal alignment files. Simply that ;)

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
