#!/usr/bin/perl -w
#
# Author: Jue Ruan <ruanjue@gmail.com>
#
use strict;

my $ref_file = shift or die("Usage: $0 <ref_file> <alignments>\n");

my %refs = ();

my $tag = undef;
my $seq = undef;
open(IN, $ref_file) or die;
while(<IN>){
	if(/^>(\S+)/){
		$refs{$tag} = $seq if(defined $tag and defined $seq);
		$tag = $1;
		$seq = '';
	} else {
		chomp;
		$seq .= $_;
	}
}
$refs{$tag} = $seq if(defined $tag and defined $seq);
close IN;

my %alns = ();

while(<>){
	my @ts = split;
	my $qry = $ts[0];
	my $obj = $ts[5];
	my $dir = $ts[6];
	my $off = $ts[8];
	my $qseq = $ts[16];
	my $tseq = $ts[17];
	push(@alns, [$qry, $obj, $dir, $off, $qseq, $tseq]);
}

foreach my $ref (sort keys %alns){
}
