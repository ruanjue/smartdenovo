#!/usr/bin/perl -w
#
# Author: Jue Ruan
#
use strict;

my @seqs = ();

while(<>){
	my @ts = split;
	push(@seqs, [$ts[0], uc $ts[1]]);
}

die("No sequences") if(@seqs < 1);

my $N = @seqs;
my $M = length $seqs[0][1];

my @weights = ();
for(my $i=0;$i<$M;$i++){
	my %bases = ();
	for(my $j=0;$j<$N;$j++){
		my $c = substr($seqs[$j][1], $i, 1);
		next if($c eq '-');
		$bases{$c} ++;
	}
	my $sum = 0;
	my $max = 0;
	foreach my $c (keys %bases){
		$sum += $bases{$c};
		$max = $bases{$c} if($bases{$c} > $max);
	}
	my $weight = 1.0;
	if($sum >= 4){
		$weight = ($sum - $max) / $max * 2.0;
	}
	push(@weights, $weight);
}

my @mats = ();
# init matches matrix
for(my $i=1;$i<$N;$i++){
	for(my $j=0;$j<$i;$j++){
		my $mat = 0;
		for(my $m=0;$m<$M;$m++){
			my $a = substr($seqs[$i][1], $m, 1);
			my $b = substr($seqs[$j][1], $m, 1);
			$mat += $weights[$m] if($a eq $b and $a ne '-');
		}
		push(@mats, $mat);
	}
}

my @tree = ();
#init nodes
for(my $i=0;$i<$N;$i++){
	push(@tree, [[], -1, 0]);
}

#clustering
for(my $iter=0;$iter<$N-1;$iter++){
	my $max = 0;
	my $mi = 0;
	my $mj = 0;
	for(my $i=0;$i<@tree;$i++){
		next if($tree[$i][1] != -1);
		for(my $j=$i+1;$j<@tree;$j++){
			next if($tree[$j][1] != -1);
			my $mat = get_upgma_matrix($i, $j);
			if($mat > $max){ $max = $mat; $mi = $i; $mj = $j; }
		}
	}
	$tree[$mi][1] = @tree;
	$tree[$mi][2] = sprintf("%0.2f", $max / 2);
	$tree[$mj][1] = @tree;
	$tree[$mj][2] = sprintf("%0.2f", $max / 2);
	push(@tree, [[$mi, $mj], -1, 0]);
	#update matches matrix
	for(my $i=0;$i<@tree-1;$i++){
		if($tree[$i][1] == -1){
			push(@mats, (get_upgma_matrix($i, $mi) + get_upgma_matrix($i, $mj)) / 2);
		} else {
			push(@mats, 0);
		}
	}
}

#print results

print_tree(scalar(@tree) - 1, 0);
print "\n"; print STDERR "\n";

1;

sub get_upgma_matrix {
	my ($i, $j) = @_;
	my $idx = ($i < $j)? ($j * ($j - 1)) / 2 + $i : ($i * ($i - 1)) / 2 + $j;
	return $mats[$idx];
}

sub print_tree {
	my $idx = shift;
	my $lv  = shift;
	my $node = $tree[$idx];
	if(@{$node->[0]}){
		print " " x $lv;
		print "(\n"; print STDERR "(";
		print_tree($node->[0][0], $lv + 1);
		print ",\n"; print STDERR ",";
		print_tree($node->[0][1], $lv + 1);
		print "\n";
		print " " x $lv;
		print "):", $node->[2]; print STDERR "):", $node->[2];
	} else {
		#print " " x $lv;
		print $seqs[$idx][0], ":", $seqs[$idx][1], ":", $node->[2];
		print STDERR $seqs[$idx][0], ":", $node->[2];
	}
}

