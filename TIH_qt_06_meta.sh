#!/usr/bin/env sh
gwas_chb=$(realpath $1)
gwas_dbds=$(realpath $2)
outdir=$(realpath $3)


module load metal/20180828


mkdir tmp/$PBS_JOBID
cd tmp/$PBS_JOBID

## SNP,CHR,BP,BETA,SE,P,N,ALLELE0,ALLELE1,A1FREQ,INFO,p_adj se_adj
echo "SNP EA OA EAF N BETA P p_adj se_adj" >input_1.txt
zcat $gwas_chb|awk 'NR>1 {print $2":"$3,$9,$8,$10,$7,$4,$6,$12,$13}' >>input_1.txt

echo "SNP EA OA EAF N BETA P p_adj se_adj" >input_2.txt
zcat $gwas_dbds|awk 'NR>1 {print $2":"$3,$9,$8,$10,$7,$4,$6,$12,$13}' >>input_2.txt


#start analysis
metal <<EOM
#input file description
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize

##=DESCRIBE AND PROCESS THE INPUT FILE ===

MARKER SNP
ALLELE EA OA
EFFECT BETA
STDERR se_adj
PVALUE p_adj
LABEL TotalSampleSize as N
FREQ EAF

PROCESS input_1.txt
PROCESS input_2.txt


##analysis
ANALYZE HETEROGENEITY
QUIT
EOM

gzip  METAANALYSIS1.TBL 
mv METAANALYSIS1.TBL.gz  $outdir/ 
mv METAANALYSIS1.TBL.info $outdir/

