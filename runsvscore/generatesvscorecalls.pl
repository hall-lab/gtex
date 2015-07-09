#!/usr/bin/perl -w

opendir(CWD, $ARGV[0]) || die "Could not open cwd: $!";
@files = grep {/^split\d\d\.vcf$/} readdir CWD;

foreach (sort @files) {
  ($prefix) = /(.*)\.vcf/;
  print "bomb -J $prefix -m 4 -e $prefix.log -o $prefix.log -q long \"/gscmnt/gc2719/halllab/users/lganel/SVScore/svscore.pl -g gencode.v19.genes.bed -e gencode.v19.exons.bed $_ > $prefix.out.vcf\"\n";
}
