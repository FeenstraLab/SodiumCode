#!/usr/bin/env sh

module load gcc/10.2.0
module load plink2/1.90beta3

## initial filter for hm3 variants and merge all chromosomes
bed_dir=$(realpath $1)
snp_dir=$(realpath $2)
rsid_dir=$(realpath $3)
sample=$(realpath $4)
output=$(realpath $5)
chr=$(echo $6)


mkdir tmp/$PBS_JOBID/
## keep id
awk 'NR>1 {print $1,$1}' $sample >tmp/$PBS_JOBID/samples.txt


plink --bfile ${bed_dir}  \
		--extract ${snp_dir}     \
		--keep tmp/$PBS_JOBID/samples.txt \
		--make-bed \
		--out tmp/$PBS_JOBID/chr${chr}_snp

plink --bfile tmp/$PBS_JOBID/chr${chr}_snp \
		--update-name ${rsid_dir} \
		--make-bed \
		--out $output

