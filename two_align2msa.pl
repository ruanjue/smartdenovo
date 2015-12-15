#!/usr/bin/perl -w
#
# Author: Jue Ruan <ruanjue@gmail.com>
#
use strict;

my @alns = ();

while(<>){
	my @ts = split;
	my $qry = $ts[0];
	my $obj = $ts[5];
	my $off = $ts[8];
	my $qseq = $ts[16];
	my $tseq = $ts[17];
	push(@alns, [$qry, $obj, $off, $qseq, $tseq]);
}

@alns = sort {$a->[2] <=> $b->[2]} @alns;

die("None alignment") unless(@alns);

my $ref_off = $alns[0][2];
my $ref_seq = '';
