#!/usr/bin/perl -w

opendir(CWD, $ARGV[0]) || die "Could not open cwd: $!";
@files = grep {/^split\d\d\.vcf$/} readdir CWD;

print "/gscmnt/gc2719/halllab/users/cchiang/src/SVScore/svscore.pl -s -g gencode.v19.genes.bed -m 4 -e gencode.v19.exons.bed -n 4 split00.vcf\n";

foreach (sort @files) {
  ($prefix) = /(.*)\.vcf/;
  print "bomb -J $prefix -m 4 -e $prefix.log -o $prefix.log -q long \"/gscmnt/gc2719/halllab/users/cchiang/src/SVScore/svscore.pl -g gencode.v19.genes.bed -m 5 -e gencode.v19.exons.bed -n 4 $_ > $prefix.out.vcf\"\n";
}
