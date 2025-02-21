#!/usr/bin/env sh

gwas_1=$(realpath $1)
gwas_2=$(realpath $2)
gwas_3=$(realpath $3)
out1=$(realpath $4)
out2=$(realpath $5)
out3=$(realpath $6)

module load anaconda2/4.4.0
module load ldsc/20200725 




mkdir tmp/$PBS_JOBID



#reformat
zcat $gwas_1|awk 'NR==1 || ($11>=0.8 && $10>=0.01 && $10<=0.99) {print $14,$9,$8,$4,$6,$7}' >tmp/$PBS_JOBID/gwas_1.txt
# SNP, A1, A2, beta, P, N
zcat $gwas_2|awk 'NR==1 || ($11>=0.8 && $10>=0.01 && $10<=0.99) {print $14,$9,$8,$4,$6,$7}' >tmp/$PBS_JOBID/gwas_2.txt


#SNP   CHR  POS_hg38 EA OA    EAF    BETA       SE       P     N  p_adj      se_adj chr
zcat $gwas_3|awk '{print $1,$4,$5,$7,$9,$10}' >tmp/$PBS_JOBID/gwas_3.txt

/tools/ldsc/20200725/munge_sumstats.py --sumstats tmp/$PBS_JOBID/gwas_1.txt --out tmp/$PBS_JOBID/gwas_1 --chunksize 500000  --snp rsid --a1 ALLELE1 --a2 ALLELE0 --p P --signed-sumstats BETA,0 --N-col N --merge-alleles /data/preprocessed/genetics/supporting_data/ldscores/w_hm3.snplist
 
/tools/ldsc/20200725/munge_sumstats.py --sumstats tmp/$PBS_JOBID/gwas_2.txt --out tmp/$PBS_JOBID/gwas_2 --chunksize 500000  --snp rsid --a1 ALLELE1 --a2 ALLELE0 --p P --signed-sumstats BETA,0 --N-col N --merge-alleles /data/preprocessed/genetics/supporting_data/ldscores/w_hm3.snplist

/tools/ldsc/20200725/munge_sumstats.py --sumstats tmp/$PBS_JOBID/gwas_3.txt --out tmp/$PBS_JOBID/gwas_3 --chunksize 500000  --snp SNP --a1 EA --a2 OA --p P --signed-sumstats BETA,0 --N-col N --merge-alleles /data/preprocessed/genetics/supporting_data/ldscores/w_hm3.snplist


## genetic correlation
/tools/ldsc/20200725/ldsc.py --rg tmp/$PBS_JOBID/gwas_1.sumstats.gz,tmp/$PBS_JOBID/gwas_2.sumstats.gz --ref-ld-chr  /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --w-ld-chr   /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --out $out1
/tools/ldsc/20200725/ldsc.py --rg tmp/$PBS_JOBID/gwas_1.sumstats.gz,tmp/$PBS_JOBID/gwas_3.sumstats.gz --ref-ld-chr  /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --w-ld-chr   /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --out $out2
/tools/ldsc/20200725/ldsc.py --rg tmp/$PBS_JOBID/gwas_2.sumstats.gz,tmp/$PBS_JOBID/gwas_3.sumstats.gz --ref-ld-chr  /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --w-ld-chr   /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --out $out3