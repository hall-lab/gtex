#!/bin/bash

# ensure that at least 10 samples are not hom ref and at least 10 samples are not hom alt

vawk --header '{ delete ALLELES; NS=0; MS=0; for (i=10;i<=NF;++i) { if (I$SVTYPE!="CNV") { if ($i!~"^0/0" && $i!~"^\\./\\.") NS+=1; if ($i!~"^1/1" && $\
                i!~"^\\./\\.") MS+=1; } else { split($i,GT,":"); ALLELES[GT[2]]+=1; } } if (I$SVTYPE=="CNV") { CN10=0; for (A in ALLELES) { if (ALLELES[A]>=10) CN10 += 1; } \
                } if ((I$SVTYPE=="CNV" && CN10>=2) || (NS>=10 && MS>=10)) print $0 }'