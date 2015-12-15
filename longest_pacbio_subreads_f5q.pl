#!/usr/bin/perl -w
#
# Author: Jue Ruan
#
use strict;

my $tag = '';
my @seqs = (['', '', '']);
my $min_len = $ENV{MIN_PB_LENGTH} || 1000;

my $line = 0;
while(<>){
	chomp;
	if(($line % 4) == 0){
		/\@(\S+)/;
		my $name = $1;
		my $subr = substr($_, length($1) + 1);
		if($name=~/^(.+?)(\/\d+_\d+)$/){
			$name = $1;
			$subr = $2 . $subr;
		}
		if($name eq $tag){
			push(@seqs, [$subr, '', '']);
		} else {
			&print_longest_seq;
			$tag = $name;
			@seqs = ([$subr, '', '']);
		}
	} elsif(($line % 4) == 1){
		$seqs[-1][1] = $_;
	} elsif(($line % 4) == 3){
		$seqs[-1][2] = $_;
	}
	$line ++;
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
	return unless($max and $max >= $min_len);
	print "\@$tag$seqs[$idx][0] len=$max\n$seqs[$idx][1]\n+\n$seqs[$idx][2]\n";
}
