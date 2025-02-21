
#!/usr/bin/env sh

module load gcc/10.2.0
module load plink2/1.90beta3
module load intel/perflibs/64 
module load regenie/3.3

bgen_file=$(realpath $1)
phe_file=$(realpath $2)
pred_file=$(realpath $3)
output_file=$(realpath $4)
CHR=$(echo $5)
sample_file=$(realpath $6)




LC_ALL=C \
LANGUAGE= \
regenie \
  --step 2 \
  --bgen ${bgen_file} \
  --sample ${sample_file} \
  --phenoFile ${phe_file} \
  --remove ${phe_file}.exclusion \
  --bsize 1000 \
  --qt \
  --chr ${CHR} \
  --pred ${pred_file} \
  --gz \
  --out ${output_file}
