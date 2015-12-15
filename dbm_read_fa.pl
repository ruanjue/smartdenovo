#!/usr/bin/perl -w
#
#Author: Ruan Jue
#
use strict;
use DB_File;

my $dbf = shift or die("Usage: $0 <dbm_file> [read_name1 ...]\n");
if($dbf!~/\.dbm$/){
	$dbf .= ".dbm" if(-e "$dbf.dbm");
}

my @names = @ARGV;

if(@names == 0){
	while(<>){
		chomp;
		push(@names, $_);
	}
}

my @tags = ();
foreach my $tag (@names){
	if($tag=~/^(.+?)(\[([+-])(:(-?\d+),(-?\d+))?\])$/){
		push(@tags, [$1, $3 eq '+'? 1:2, (defined $5)? $5:1, (defined $6)? $6:-1, 1]);
	} else {
		push(@tags, [$tag, 1, 1, -1, 0]);
	}
}

my %seqs;

tie %seqs, 'DB_File', $dbf, O_RDONLY or die "Cannot open $dbf: $!";

foreach my $tag (@tags){
	if(exists $seqs{$tag->[0]}){
		my $seq = $seqs{$tag->[0]};
		$tag->[3] = length($seq) if($tag->[3] < 1);
		if($tag->[4]){
			print ">", join("_", $tag->[0], $tag->[1] == 1? "F":"R", $tag->[2], $tag->[3]), "\n";
		} else {
			print ">$tag->[0]\n";
		}
		if($tag->[2] < $tag->[3]){
			my $ss = substr($seq, $tag->[2] - 1, $tag->[3] - $tag->[2] + 1);
			if($tag->[1] == 2){
				$ss =~tr/ACGTacgt/TGCAtgca/;
				$ss = reverse $ss;
			}
			while($ss=~/(.{1,100})/g){ print "$1\n"; }
		}
	} else {
		warn("'$tag->[0]' was not found\n");
	}
}

untie %seqs;

1;
