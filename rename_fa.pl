#!/usr/bin/perl -w
#
# Author: Ruan Jue
#
use strict;
use Getopt::Std;

my $prefix = $ENV{'PARAM_RENAME_FA_PREFIX'} || 'S';
my $suffix = $ENV{'PARAM_RENAME_FA_SUFFIX'} || '';

our ($opt_p, $opt_s, $opt_h, $opt_f);

getopts("hp:s:f:");
die("Usage: $0 [-p name_prefix] [-s name_suffix] [-f trans_file] <fasta_file>\n") if($opt_h);
$prefix = $opt_p if(defined $opt_p);
$suffix = $opt_s if(defined $opt_s);
my %hash;
if(defined $opt_f){
	open(IN, "<", $opt_f) or die;
	%hash = ();
	while(<IN>){
		my @ts = split;
		$hash{$ts[0]} = $ts[1];
		print STDERR "$ts[0]\t$ts[1]\n";
	}
	close IN;
}

my $idx = 0;

while(<>){
	if(/^>(\S+)/){
		my $desc = substr($_, length($1));
		$idx ++;
		if(%hash){
			if(exists $hash{$1}){
				my $tag = $hash{$1};
				print ">$tag $1$desc", substr($_, length($1) + 1);
			} else {
				print;
			}
		} else {
			printf(">$prefix%012d$suffix $1$desc", $idx);
		}
	} else {
		print;
	}
}

1;
