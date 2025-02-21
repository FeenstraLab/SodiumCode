#!/usr/bin/env sh

gwas_1=$(realpath $1)
gwas_2=$(realpath $2)
out1=$(realpath $3)

module load anaconda2/4.4.0
module load ldsc/20200725 




mkdir tmp/$PBS_JOBID



#reformat to SNP, A1, A2, beta, P, N

#CHR BP chr MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Effect  StdErr  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal TotalSampleSize rsid
zcat $gwas_1|awk '{print $20,$5,$6,$11,$13,$19}' >tmp/$PBS_JOBID/gwas_1.txt


#SNP   CHR  POS_hg38 EA OA    EAF    BETA       SE       P     N  p_adj      se_adj chr
zcat $gwas_2|awk '{print $1,$4,$5,$7,$9,$10}' >tmp/$PBS_JOBID/gwas_2.txt

/tools/ldsc/20200725/munge_sumstats.py --sumstats tmp/$PBS_JOBID/gwas_1.txt --out tmp/$PBS_JOBID/gwas_1 --chunksize 500000  --snp rsid --a1 Allele1 --a2 Allele2 --p P.value --signed-sumstats Effect,0 --N-col TotalSampleSize --merge-alleles /data/preprocessed/genetics/supporting_data/ldscores/w_hm3.snplist
 
/tools/ldsc/20200725/munge_sumstats.py --sumstats tmp/$PBS_JOBID/gwas_2.txt --out tmp/$PBS_JOBID/gwas_2 --chunksize 500000  --snp SNP --a1 EA --a2 OA --p P --signed-sumstats BETA,0 --N-col N --merge-alleles /data/preprocessed/genetics/supporting_data/ldscores/w_hm3.snplist


## genetic correlation
/tools/ldsc/20200725/ldsc.py --rg tmp/$PBS_JOBID/gwas_1.sumstats.gz,tmp/$PBS_JOBID/gwas_2.sumstats.gz --ref-ld-chr  /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --w-ld-chr   /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --out $out1