#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (t=>8, k=>17, m=>0.6, s=>200, J=>0, p=>'wtasm');
getopts('t:p:k:m:s:J:', \%opts);
die (qq/Usage: smartdenovo.pl [options] <reads.fa>
Options:
  -p STR     output prefix [$opts{p}]
  -t INT     number of threads [$opts{t}]
  -k INT     k-mer length for overlapping [$opts{k}]
  -J INT     min read length [$opts{J}]
  -m FLOAT   min identity [$opts{m}]
  -s INT     min alignment score [$opts{s}]
/) if (@ARGV == 0);

$opts{pre} = gwhich("wtpre") || die;
$opts{zmo} = gwhich("wtzmo") || die;
$opts{gbo} = gwhich("wtgbo") || die;
$opts{clp} = gwhich("wtclp") || die;
$opts{lay} = gwhich("wtlay") || die;
$opts{cns} = gwhich("wtcns") || die;

my @lines = ();
push(@lines, qq/PREFIX=$opts{p}/, '');
push(@lines, qq/EXE_PRE=$opts{pre}/, qq/EXE_ZMO=$opts{zmo}/, qq/EXE_GBO=$opts{gbo}/, qq/EXE_CLP=$opts{clp}/, qq/EXE_LAY=$opts{lay}/, qq/EXE_CNS=$opts{cns}/);
push(@lines, qq/N_THREADS=$opts{t}/, "");

push(@lines, q/all:$(PREFIX).cns/, "");

push(@lines, q/$(PREFIX).pre.gz:/);
push(@lines, qq/\t\$(EXE_PRE) -L -J $opts{J} $ARGV[0] 2> \$\@.log | gzip -1 > \$\@/, "");

push(@lines, q/$(PREFIX).ovl.short:$(PREFIX).pre.gz/);
push(@lines, qq/\t\$(EXE_ZMO) -t \$(N_THREADS) -i \$(PREFIX).pre.gz -fo - -k $opts{k} -s $opts{s} -m $opts{m} 2> \$\@.log | cut -f1-16 > \$\@/, "");

push(@lines, q/$(PREFIX).gbo.short:$(PREFIX).pre.gz $(PREFIX).ovl.short/);
push(@lines, qq/\t\$(EXE_GBO) -t \$(N_THREADS) -i \$< -j \$(PREFIX).ovl.short -fo - 2> \$\@.log | cut -f1-16 > \$\@/, "");

push(@lines, q/$(PREFIX).obt:$(PREFIX).ovl.short $(PREFIX).gbo.short/);
push(@lines, qq/\t\$(EXE_CLP) -i \$< -i \$(PREFIX).gbo.short -fo \$\@ -F -d 4 2> \$\@.log/, "");

push(@lines, q/$(PREFIX).lay.utg $(PREFIX).lay:$(PREFIX).pre.gz $(PREFIX).ovl.short $(PREFIX).gbo.short $(PREFIX).obt/);
push(@lines, qq/\t\$(EXE_LAY) -i \$< -b \$(PREFIX).obt -j \$(PREFIX).ovl.short -j \$(PREFIX).gbo.short -fo \$(PREFIX).lay -s $opts{s} -m $opts{m} -R -r 1 -c 1 2> \$(PREFIX).lay.log/, "");

push(@lines, q/$(PREFIX).cns:$(PREFIX).lay/);
push(@lines, qq/\t\$(EXE_CNS) -t \$(N_THREADS) \$< > \$\@ 2> \$\@.log/, "");

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
