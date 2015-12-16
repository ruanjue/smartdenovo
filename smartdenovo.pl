#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (t=>8, k=>17, m=>0.6, s=>200, p=>'wtasm');
getopts('t:p:k:m:s:', \%opts);
die (qq/Usage: smartdenovo.pl [options] <reads.fa>
Options:
  -p STR     output prefix [$opts{p}]
  -t INT     number of threads [$opts{t}]
  -k INT     k-mer length for overlapping [$opts{k}]
  -m FLOAT   min identity [$opts{m}]
  -s INT     min alignment score [$opts{s}]
/) if (@ARGV == 0);

$opts{zmo} = gwhich("wtzmo") || die;
$opts{gbo} = gwhich("wtgbo") || die;
$opts{clp} = gwhich("wtclp") || die;
$opts{lay} = gwhich("wtlay") || die;

my $input = $ARGV[0];

my @lines = ();
push(@lines, qq/PREFIX=$opts{p}/, '');
push(@lines, qq/EXE_ZMO=$opts{zmo}/, qq/EXE_GBO=$opts{gbo}/, qq/EXE_CLP=$opts{clp}/, qq/EXE_LAY=$opts{lay}/);
push(@lines, qq/N_THREADS=$opts{t}/, "");
push(@lines, qq/INPUT=$input/, "");

push(@lines, qq/all:\$(PREFIX).lay.utg/, "");

push(@lines, q/$(PREFIX).ovl.short:/);
push(@lines, qq/\t\$(EXE_ZMO) -t \$(N_THREADS) -i \$(INPUT) -fo - -k $opts{k} -s $opts{s} -m $opts{m} 2> \$\@.log | cut -f1-16 > \$\@/, "");

push(@lines, q/$(PREFIX).gbo.short:$(PREFIX).ovl.short/);
push(@lines, qq/\t\$(EXE_GBO) -t \$(N_THREADS) -i \$(INPUT) -j \$< -fo - 2> \$\@.log | cut -f1-16 > \$\@/, "");

push(@lines, q/$(PREFIX).obt:$(PREFIX).ovl.short $(PREFIX).gbo.short/);
push(@lines, qq/\t\$(EXE_CLP) -i \$< -i \$(PREFIX).gbo.short -fo \$\@ -F -d 4 2> \$\@.log/, "");

push(@lines, q/$(PREFIX).lay.utg:$(PREFIX).ovl.short $(PREFIX).gbo.short $(PREFIX).obt/);
push(@lines, qq/\t\$(EXE_LAY) -i \$(INPUT) -b \$(PREFIX).obt -j \$< -j \$(PREFIX).gbo.short -fo \$(PREFIX).lay -s $opts{s} -m $opts{m} -R -r 1 -c 1 2> \$(PREFIX).lay.log/, "");

print(join("\n", @lines), "\n");

sub which {
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file") && (-f "$x/$file");
	}
	return;
}

sub gwhich {
    my $progname = shift;
    my $addtional_path = shift if (@_);
    my $dirname = &dirname($0);
    my $tmp;

    chomp($dirname);
    if ($progname =~ /^\// && (-x $progname) && (-f $progname)) {
        return $progname;
    } elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
        return $tmp;
    } elsif (defined($dirname) && (-x "$dirname/$progname") && (-f "$dirname/$progname")) {
        return "$dirname/$progname";
    } elsif ((-x "./$progname") && (-f "./$progname")) {
        return "./$progname";
    } elsif (($tmp = &which($progname))) {
        return $tmp;
    } else {
        return;
    }
}

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
