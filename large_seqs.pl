#!/usr/bin/perl -w
#
#Author: Ruan Jue<ruanjue@genomics.org.cn>

use warnings;
use strict;

my $rev = $ENV{REVERSE} || 0;
my $seq_file = shift or die("Usage: [REVERSE=0|1] $0 <fasta_file> [seq_len:100]\n");
my $seq_len  = shift || 100;

my $doc = '';
my $seq = '';
my $len = 0;

my $in;
if($seq_file eq'-'){
	$in = \*STDIN;
} else {
	open($in, $seq_file) or die($!);
}

while(<$in>){
	if(/^>/){
		my $line = $_;
		if($doc and ($rev? $len <= $seq_len : $len >= $seq_len)){
			print $doc;
			print $seq;
		}
		$doc = $line;
		$seq = '';
		$len = 0;
	} else {
		$len += length($_) - 1;
		$seq .= $_;
	}
}
if($doc and ($rev? $len <= $seq_len : $len >= $seq_len)){
	print $doc;
	print $seq;
}
close $in if($seq_file ne '-');
