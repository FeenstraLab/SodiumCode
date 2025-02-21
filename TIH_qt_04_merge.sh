#!/usr/bin/env sh

inputdir=$(realpath $1)
outputdir=$(realpath $2)

echo "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA P" > ${outputdir}/gwas.txt

for CHR in $(seq 1 22) X
do
                #CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
                zcat ${inputdir}/step2_chr${CHR}_phe.regenie.gz|awk 'NR>1 {print $0,10^(-$13)}' >> ${outputdir}/gwas.txt
done

wc -l ${outputdir}/gwas.txt
gzip ${outputdir}/gwas.txt
