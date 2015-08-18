#!/usr/bin/perl -w

opendir(CWD, $ARGV[0]) || die "Could not open cwd: $!";
@files = grep {/^split\d\d\.vcf$/} readdir CWD;

print "/gscmnt/gc2719/halllab/users/lganel/SVScore/svscore.pl -s -g /gscmnt/gc2719/halllab/users/lganel/SVScore/gencode.v19.genes.bed -e /gscmnt/gc2719/halllab/users/lganel/SVScore/gencode.v19.exons.bed -n 4 -c /gscmnt/gc2719/halllab/users/lganel/SVScore/whole_genome_SNVs.tsv.gz -o max split00.vcf\n";

foreach (sort @files) {
  ($prefix) = /(.*)\.vcf/;
  print "bomb -J $prefix -m 4 -e $prefix.log -o $prefix.log -q hall-lab \"/gscmnt/gc2719/halllab/users/lganel/SVScore/svscore.pl -g /gscmnt/gc2719/halllab/users/lganel/SVScore/gencode.v19.genes.bed -e /gscmnt/gc2719/halllab/users/lganel/SVScore/gencode.v19.exons.bed -n 4 -c /gscmnt/gc2719/halllab/users/lganel/SVScore/whole_genome_SNVs.tsv.gz -o max $_ > $prefix.out.vcf\"\n";
}


# /gscmnt/gc2719/halllab/users/cchiang/src/SVScore/svscore.pl -d -g gencode.v19.genes.bed -e gencode.v19.exons.bed -n 4  -c /gscuser/lganel/lganel/SVScore/whole_genome_SNVs.tsv.gz test.vcf
