#!/usr/bin/env sh

module load gcc/10.2.0
module load plink2/1.90beta3
module load intel/perflibs/64 
module load regenie/2.2.4

geno_file=$(realpath $1)
phe_file=$(realpath $2)
cov_file=$(realpath $3)
output_file=$(realpath $4)





LC_ALL=C \
LANGUAGE= \
regenie \
--step 1 \
--bed ${geno_file} \
--covarFile ${cov_file} \
--phenoFile ${phe_file} \
--covarColList thiazide_age,PC{1:5} \
--catCovarList sex,dbds \
--bsize 1000 \
--bt \
--lowmem \
--lowmem-prefix ${output_file} \
--gz \
--out ${output_file}
