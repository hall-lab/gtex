bomb -J split00 -m 4 -e split00.log -o split00.log -q long "/gscmnt/gc2719/halllab/users/lganel/SVScore/svscore.pl -g gencode.v19.genes.bed -e gencode.v19.exons.bed split00.vcf > split00.out.vcf"
bomb -J split01 -m 4 -e split01.log -o split01.log -q long "/gscmnt/gc2719/halllab/users/lganel/SVScore/svscore.pl -g gencode.v19.genes.bed -e gencode.v19.exons.bed split01.vcf > split01.out.vcf"
