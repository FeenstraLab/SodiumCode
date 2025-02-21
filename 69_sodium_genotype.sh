#!/usr/bin/env sh

module load gcc/10.2.0
module load plink2/1.90beta3

## initial filter for hm3 variants and merge all chromosomes
bed_dir=$(realpath $1)
snp_dir=$(realpath $2)
rsid_dir=$(realpath $3)
lab=$(realpath $4)
output=$(realpath $5)
chr=$(echo $6)


mkdir tmp/$PBS_JOBID/
## keep id
zgrep  "NPU03429" $lab | awk '{print $1}'|sort |uniq| awk '{print $1,$1}'  >tmp/$PBS_JOBID/samples.txt

# for chr in $(seq 1 22)
# do
plink --bfile ${bed_dir}/20210503_chr${chr}_CVD_FINAL  \
		--extract ${snp_dir}/chr${chr}_hm3_SNP.txt     \
		--keep tmp/$PBS_JOBID/samples.txt \
		--make-bed \
		--out tmp/$PBS_JOBID/chr${chr}_snp \
		--threads 8

plink --bfile tmp/$PBS_JOBID/chr${chr}_snp \
		--update-name ${rsid_dir}/chr${chr}_hm3_rsid.tsv \
		--make-bed \
		--out $output/chr${chr} \
		--threads 8
# done

# echo "tmp/$PBS_JOBID/chr2
# tmp/$PBS_JOBID/chr3
# tmp/$PBS_JOBID/chr4
# tmp/$PBS_JOBID/chr5
# tmp/$PBS_JOBID/chr6
# tmp/$PBS_JOBID/chr7
# tmp/$PBS_JOBID/chr8
# tmp/$PBS_JOBID/chr9
# tmp/$PBS_JOBID/chr10
# tmp/$PBS_JOBID/chr11
# tmp/$PBS_JOBID/chr12
# tmp/$PBS_JOBID/chr13
# tmp/$PBS_JOBID/chr14
# tmp/$PBS_JOBID/chr15
# tmp/$PBS_JOBID/chr16
# tmp/$PBS_JOBID/chr17
# tmp/$PBS_JOBID/chr18
# tmp/$PBS_JOBID/chr19
# tmp/$PBS_JOBID/chr20
# tmp/$PBS_JOBID/chr21
# tmp/$PBS_JOBID/chr22" >tmp/$PBS_JOBID/chr2_chr22.list

# plink --bfile tmp/$PBS_JOBID/chr1 --merge-list tmp/$PBS_JOBID/chr2_chr22.list --make-bed --out $output --threads 8
