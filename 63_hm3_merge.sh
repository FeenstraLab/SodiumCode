#!/usr/bin/env sh

module load gcc/10.2.0
module load plink2/1.90beta3

## initial filter for hm3 variants and merge all chromosomes
input_dir=$(realpath $1)
output=$(realpath $2)



mkdir tmp/$PBS_JOBID/


echo "${input_dir}/chr2
${input_dir}/chr3
${input_dir}/chr4
${input_dir}/chr5
${input_dir}/chr6
${input_dir}/chr7
${input_dir}/chr8
${input_dir}/chr9
${input_dir}/chr10
${input_dir}/chr11
${input_dir}/chr12
${input_dir}/chr13
${input_dir}/chr14
${input_dir}/chr15
${input_dir}/chr16
${input_dir}/chr17
${input_dir}/chr18
${input_dir}/chr19
${input_dir}/chr20
${input_dir}/chr21
${input_dir}/chr22" >tmp/$PBS_JOBID/chr2_chr22.list

plink --bfile ${input_dir}/chr1 --merge-list tmp/$PBS_JOBID/chr2_chr22.list --make-bed --out $output


## allele frequency
plink --bfile $output --freq --out $output