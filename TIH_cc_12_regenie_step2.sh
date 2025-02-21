
#!/usr/bin/env sh

module load gcc/12.2.0
module load plink2/1.90beta3
module load intel/perflibs/64 
module load regenie/2.2.4

bgen_file=$(realpath $1)
phe_file=$(realpath $2)
cov_file=$(realpath $3)
pred_file=$(realpath $4)
output_file=$(realpath $5)
CHR=$(echo $6)
sample_file=$(realpath $7)




LC_ALL=C \
LANGUAGE= \
/services/tools/regenie/3.3/regenie \
  --step 2 \
  --bgen ${bgen_file} \
  --sample ${sample_file} \
  --covarFile ${cov_file} \
  --phenoFile ${phe_file} \
  --covarColList thiazide_age,PC{1:5} \
  --catCovarList sex,dbds \
  --bsize 1000 \
  --bt \
  --chr ${CHR} \
  --pred ${pred_file} \
  --gz \
  --out ${output_file}
