#!/usr/bin/env sh


gwas_1=$(realpath $1)
out1=$(realpath $2)

## meta-analysis heritability check

## extract hm3 variants and get rsid
#	module load gcc/9.3.0 perl/5.30.2 intel/redist/2019_update2 intel/perflibs/2019_update2 openssl/1.1.1b  apache-arrow/11.0.0_CPP
#	module load R/4.2.0

#	Rscript scripts/d20230804/h2_check.R



#### format
module load anaconda2/4.4.0
module load ldsc/20200725 


#gwas_1="results/d20230804/meta_1w/excl3_lm_le50_weights/gwas_hm3.txt.gz"
#out1="results/d20230804/meta_1w/excl3_lm_le50_weights/h2"

mkdir tmp/$PBS_JOBID

	# format to SNP, A1, A2, beta, P, N
	zcat $gwas_1|awk '{print $18,$3,$4,$9,$11,$17}' >tmp/$PBS_JOBID/gwas_1.txt

	/tools/ldsc/20200725/munge_sumstats.py --sumstats tmp/$PBS_JOBID/gwas_1.txt --out tmp/$PBS_JOBID/gwas_1 --chunksize 500000  --snp rsid --a1 Allele1 --a2 Allele2 --p P --signed-sumstats Effect,0 --N-col TotalSampleSize --merge-alleles /data/preprocessed/genetics/supporting_data/ldscores/w_hm3.snplist

	#h2
	/tools/ldsc/20200725/ldsc.py --h2 tmp/$PBS_JOBID/gwas_1.sumstats.gz --ref-ld-chr  /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --w-ld-chr   /data/preprocessed/genetics/supporting_data/ldsc/eur_w_ld_chr_2021/ --out $out1
