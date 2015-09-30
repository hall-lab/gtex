#!/bin/bash -e

USAGE='\tusage: '$0' <VCF> <outliers> <distance>'

if [ $# -ne 3 ]
then
    echo -e "\n$USAGE\n"
    exit
fi

VCF=$1
OUTLIERS=$2
DISTANCE=$3

# only look for private and doubleton
NSAMP_MAX=2

while read line
do
    # echo $line

    SAMPLE=`echo "$line" | awk '{ print $2 }'`
    GENE=`echo "$line" | awk '{ print $3 }'`
    TISSUE=`echo "$line" | awk '{ print $4 }'`
    Z=`echo "$line" | awk '{ print $5 }'`
    ID=`echo "$line" | awk '{ print $1 }'`

    GENE_COORD=(`zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.genes.bed.gz | awk -v GENE=$GENE '$4==GENE' | sort -k1,1V -k2,2n | bedtools merge | head -n 1 | awk -v DISTANCE=$DISTANCE '{ print $1,$2,$3 }'`)
    REGION=`echo -e "${GENE_COORD[0]}\t${GENE_COORD[1]}\t${GENE_COORD[2]}" | awk -v DISTANCE=$DISTANCE '{ print $1":"$2-DISTANCE"-"$3+DISTANCE }'`

    VARS=`tabix -h $VCF $REGION \
	| vawk --header -v REGION=$REGION '{ split(REGION,r1,":"); split(r1[2],r2,"-"); CHROM=r1[1]; START=r2[1]; STOP=r2[2]; if (I$SVTYPE=="INV") { if ($1==CHROM && (($2>START && $2<STOP) || (I$END>START && I$END<STOP))) print } else print }' \
	| vawk --header "{ if (S\\$$SAMPLE ~ \"^0/1\" || S\\$$SAMPLE ~ \"^1/1\") print }" \
	| vawk -v S=$SAMPLE -v GENE=$GENE -v TISSUE=$TISSUE -v Z=$Z -v ID=$ID -v MAX_NSAMP=$NSAMP_MAX -v GENE_CHROM=${GENE_COORD[0]} -v GENE_START=${GENE_COORD[1]} -v GENE_END=${GENE_COORD[2]} '{ if (I$NSAMP<=MAX_NSAMP) print $1,$2,$3,$4,$5,$6,$7,$8,ID,S,GENE,TISSUE,Z,GENE_CHROM,GENE_START,GENE_END }'`

    MATCH=0
    if [[ ! -z "$VARS" ]]
    then
	MATCH=1
    fi
    if [[ "$MATCH" -eq 1 ]]
    then
	echo "$VARS"
    fi
done < <(zless $OUTLIERS)
