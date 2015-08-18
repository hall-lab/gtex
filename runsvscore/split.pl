#!/usr/bin/perl -w

##sort by ID first!!
## grep -v '^#' gtex_merged.sv.gt.cn.class.af_gt_0.ann.vcf | sort -k3,3 | cat <(grep '^#' gtex_merged.sv.gt.cn.class.af_gt_0.ann.vcf) - | perl split.pl

$linenum = 1;
$filenum = 0;
$linesperfile = 800; ##CHANGEME

open(OUT,"> split00.vcf") || die "Could not open output file 00: $!";

while (<>) {
  if (/^#/) {
    push @header,$_;
    next;
  }

  if ($linenum == 1) {
    foreach (@header) {
      print OUT;
    }
  }

  if ($linenum >= $linesperfile && /SVTYPE=BND/) {
    @split = (split(/\s+/));
    @lastsplit = split(/\s+/,$lines[$#lines]);
    ($currentid) = ($split[2] =~ /^(\d+)_[12]$/);
    ($lastid) = ($lastsplit[2] =~ /^(\d+)(?:_[12])?$/);
    $split = ($currentid == $lastid);
  } else {
    $split = ($linenum >= $linesperfile);
  }
  push @lines, $_;
  
  if ($split) {
    $split = 0;
    $linenum = 0;
    %coords = map {$_ => [(split(/\s+/,$_))[0..1]]} @lines;
    foreach (sort {$coords{$a}->[0] cmp $coords{$b}->[0] || $coords{$a}->[1] cmp $coords{$b}->[1]} @lines) {
      print OUT;
    }
    @lines = ();
    close OUT;
    $filenum++;
    $filestring = ($filenum < 10 ? "0$filenum" : $filenum);
    open(OUT, "> split$filestring.vcf") || die "Could not open output file $filestring: $!";
  }

  $linenum++;

}

%coords = map {$_ => [(split(/\s+/,$_))[0..1]]} @lines;
foreach (sort {$coords{$a}->[0] cmp $coords{$b}->[0] || $coords{$a}->[1] cmp $coords{$b}->[1]} @lines) {
  print OUT;
}
