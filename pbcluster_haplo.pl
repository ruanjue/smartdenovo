#!/usr/bin/perl -w
#
# AUthor: Jue Ruan
#
use strict;

my $core = shift or die("Usage: $0 <pb_name> <msa>\n");

my @seqs = ();
my $ref = undef;

while(<>){
	my @ts = split;
	$ts[1] = uc $ts[1];
	push(@seqs, [$ts[0], $ts[1], 0]);
	if($ts[0] eq $core){
		$ref = $ts[1];
	}
}

die("No sequences") if(@seqs < 1);
die("Cannot find $core") if(not defined $ref);

my $N = @seqs;
my $M = length $ref;
my @weights = ();
for(my $i=0;$i<$M;$i++){ push(@weights, 1); }

for(my $iter=0;$iter<4;$iter++){

	for(my $i=0;$i<$N;$i++){
		my $mat = 0;
		for(my $j=0;$j<$M;$j++){
			my $a = substr($ref, $j, 1);
			my $b = substr($seqs[$i][1], $j, 1);
			$mat += $weights[$j] if($a eq $b and $a ne '-');
		}
		$seqs[$i][2] = $mat;
	}

	@seqs = sort {$b->[2] <=> $a->[2]} @seqs;

	my $par = 1;
	my $lst = -1;
	while(1){
		my $cnt = int($N / $par);
		$par ++;
		last if($cnt < 3);
		next if($cnt == $lst);
		$lst = $cnt;
		call_cns($cnt);
		#last if($par > 8);
	}

	print join("\n", map {join("\t", @$_)} @seqs), "\n";

}

1;

sub call_cns {
	my $CNT = shift;
	my @cns = ();
	for(my $i=0;$i<$M;$i++){
		push(@cns, '-'),next if(substr($ref, $i, 1) eq '-');
		my @bases = ();
		for(my $j=0;$j<$N;$j++){
			my $c = substr($seqs[$j][1], $i, 1);
			push(@bases, [$c, $N - $j]) if($c ne '-');
		}
		my %hash = ();
		my $cnt = $CNT;
		$cnt = @bases if($cnt > @bases);
		for(my $j=0;$j<$cnt;$j++){
			my $c = $bases[$j];
			next if($c->[0] eq '-');
			$hash{$c->[0]} += $c->[1];
		}
		my $max = [1, '-'];
		foreach my $c (sort keys %hash){
			#print "$c:$hash{$c}\t";
			$max = [$hash{$c}, $c] if($hash{$c} > $max->[0]);
		}
		#print "\n";
		#print "$i\t$max->[1]\t", join('', map{$_->[0]} @bases), "\n";
		push(@cns, $max->[1]);
	}

	$ref = join('', @cns);
	print "REF[$CNT]\t$ref\n";
}

