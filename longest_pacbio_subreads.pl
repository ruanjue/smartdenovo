#!/usr/bin/perl -w
#
# Author: Jue Ruan
#
use strict;

my $tag = '';
my @seqs = (['', '']);

while(<>){
	if(/^>(\S+)/){
		my $name = $1;
		my $subr = '';
		if($name=~/^(.+?)(\/\d+_\d+)$/){
			$name = $1;
			$subr = $2;
		}
		if($name eq $tag){
			push(@seqs, [$subr, '']);
		} else {
			&print_longest_seq;
			$tag = $name;
			@seqs = ([$subr, '']);
		}
	} else {
		chomp;
		$seqs[-1][1] .= $_;
	}
}
&print_longest_seq;

1;

sub print_longest_seq {
	my $idx = 0;
	my $max = 0;
	for(my $i=0;$i<@seqs;$i++){
		if(length($seqs[$i][1]) > $max){
			$idx = $i; $max = length($seqs[$i][1]);
		}
	}
	return unless($max);
	print ">$tag$seqs[$idx][0] len=$max\n$seqs[$idx][1]\n";
}
