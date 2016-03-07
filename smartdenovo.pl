#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (t=>8, e=>'dmo', k=>16, J=>5000, p=>'wtasm', c=>0);
getopts('t:p:e:k:m:s:J:c:', \%opts);
$opts{e} = lc $opts{e};
die (qq/Usage: smartdenovo.pl [options] <reads.fa>
Options:
  -p STR     output prefix [$opts{p}]
  -e STR     engine of overlaper, compressed kmer overlapper(zmo), dot matrix overlapper(dmo), [$opts{e}]
  -t INT     number of threads [$opts{t}]
  -k INT     k-mer length for overlapping [$opts{k}]
             for large genome as human, please use 17
  -J INT     min read length [$opts{J}]
  -c INT     generate consensus, [$opts{c}]
/) if (@ARGV == 0);

die("Engine choice: dmo, zmo\n") if($opts{e} ne 'dmo' and $opts{e} ne 'zmo');

$opts{pre} = gwhich("wtpre") || die;
$opts{zmo} = gwhich("wtzmo") || die;
$opts{obt} = gwhich("wtobt") || die;
$opts{gbo} = gwhich("wtgbo") || die;
$opts{clp} = gwhich("wtclp") || die;
$opts{lay} = gwhich("wtlay") || die;
$opts{cns} = gwhich("wtcns") || die;


my @lines = ();
push(@lines, qq/PREFIX=$opts{p}/, '');
push(@lines, qq/EXE_PRE=$opts{pre}/, qq/EXE_ZMO=$opts{zmo}/, qq/EXE_OBT=$opts{obt}/, qq/EXE_GBO=$opts{gbo}/, qq/EXE_CLP=$opts{clp}/, qq/EXE_LAY=$opts{lay}/, qq/EXE_CNS=$opts{cns}/);
push(@lines, qq/N_THREADS=$opts{t}/, "");

if($opts{c}){
	push(@lines, qq/all:\$(PREFIX).$opts{e}.cns/, "");
} else {
	push(@lines, qq/all:\$(PREFIX).$opts{e}.lay/, "");
}

push(@lines, q/$(PREFIX).fa.gz:/);
push(@lines, qq/\t\$(EXE_PRE) -J $opts{J} $ARGV[0] | gzip -c -1 > \$\@/, "");

if($opts{e} eq 'dmo'){
	push(@lines, q/$(PREFIX).dmo.ovl:$(PREFIX).fa.gz/);
	push(@lines, qq/\t\$(EXE_ZMO) -t \$(N_THREADS) -i \$(PREFIX).fa.gz -fo \$\@ -k $opts{k} -z 10 -Z 16 -U -1 -m 0.1 -A 1000/, "");

	push(@lines, q/$(PREFIX).dmo.obt:$(PREFIX).fa.gz $(PREFIX).dmo.ovl/);
	push(@lines, qq/\t\$(EXE_CLP) -i \$(PREFIX).dmo.ovl -fo \$\@ -d 3 -k 300 -m 0.1 -FT/, "");

	push(@lines, q/$(PREFIX).dmo.lay:$(PREFIX).fa.gz $(PREFIX).dmo.obt $(PREFIX).dmo.ovl/);
	push(@lines, qq/\t\$(EXE_LAY) -i \$(PREFIX).fa.gz -b \$(PREFIX).dmo.obt -j \$(PREFIX).dmo.ovl -fo \$(PREFIX).dmo.lay -w 300 -s 200 -m 0.1 -r 0.95 -c 1/, "");

} else {
	push(@lines, q/$(PREFIX).zmo.ovl.short:$(PREFIX).fa.gz/);
	push(@lines, qq/\t\$(EXE_ZMO) -t \$(N_THREADS) -i \$(PREFIX).fa.gz -fo - -k $opts{k} -s 200 -m 0.6 | cut -f1-16 > \$\@/, "");

	push(@lines, q/$(PREFIX).zmo.gbo.short:$(PREFIX).fa.gz $(PREFIX).zmo.ovl.short/);
	push(@lines, qq/\t\$(EXE_GBO) -t \$(N_THREADS) -i \$< -j \$(PREFIX).zmo.ovl.short -fo - | cut -f1-16 > \$\@/, "");

	push(@lines, q/$(PREFIX).zmo.obt:$(PREFIX).zmo.ovl.short $(PREFIX).zmo.gbo.short/);
	push(@lines, qq/\t\$(EXE_CLP) -i \$< -i \$(PREFIX).zmo.gbo.short -fo \$\@ -F -d 2/, "");

	push(@lines, q/$(PREFIX).zmo.lay:$(PREFIX).fa.gz $(PREFIX).zmo.ovl.short $(PREFIX).zmo.gbo.short $(PREFIX).zmo.obt/);
	push(@lines, qq/\t\$(EXE_LAY) -i \$< -b \$(PREFIX).zmo.obt -j \$(PREFIX).zmo.ovl.short -j \$(PREFIX).zmo.gbo.short -fo \$(PREFIX).zmo.lay -s 200 -m 0.6 -R -r 1 -c 1/, "");
}

if($opts{c}){
	push(@lines, qq/\$(PREFIX).$opts{e}.cns:\$(PREFIX).$opts{e}.lay/);
	push(@lines, qq/\t\$(EXE_CNS) -t \$(N_THREADS) \$< > \$\@ 2> \$\@.log/, "");
}
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
