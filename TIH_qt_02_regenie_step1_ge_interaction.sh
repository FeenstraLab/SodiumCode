#!/usr/bin/env sh

module load gcc/12.2.0
module load plink2/1.90beta3
module load intel/perflibs/64 
module load regenie/3.3

geno_file=$(realpath $1)
phe_file=$(realpath $2)
output_file=$(realpath $3)
exclusion_list=$(realpath $4)

awk 'NR>1' $exclusion_list >${phe_file}.exclusion



LC_ALL=C \
LANGUAGE= \
regenie \
--step 1 \
--bed ${geno_file} \
--phenoFile ${phe_file} \
--remove ${phe_file}.exclusion \
--bsize 1000 \
--qt \
--lowmem \
--lowmem-prefix ${output_file} \
--gz \
--out ${output_file}


