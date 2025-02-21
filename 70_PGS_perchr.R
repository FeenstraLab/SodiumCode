library(optparse)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bigreadr)

parsed_opts <- list(
  make_option("--bed",
              action = "store",
              default = NA,
              type = 'character',
              help = 'bed file for each chromosome'),
  make_option(c("-t", "--threads"),
              action = "store",
              default = 1,
              type = 'integer',
              help = "Number of threads to use"),
  make_option("--beta_auto",
              action = "store",
              default = NA,
              type = 'character',
              help = 'input of LDpred2 auto model beta effect'),
    make_option("--out_pred",
              action = "store",
              default = NA,
              type = 'character',
              help = ' PRS output')
  
)

if (interactive()) {


} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

## input matrix
rds <- snp_readBed2(opt$bed,backingfile = tempfile(),ncores = opt$threads)


obj.bigSNP <- snp_attach(rds)
G <- obj.bigSNP$genotypes

## input auto mdoel results
beta_auto <-  readRDS(opt$beta_auto) 


map <- obj.bigSNP$map %>%
        rename(chr=chromosome,rsid=marker.ID,pos=physical.pos,a1=allele1,a0=allele2)

info_snps <- snp_match(beta_auto,map,join_by_pos = FALSE)

# info_snps <- snp_match(sumstats, map,join_by_pos = FALSE)


## impute G
print("impute G")
print(format(Sys.time(), "%X"))
G_imp <- snp_fastImputeSimple(
  G,
  ncores = opt$threads
)
print(format(Sys.time(), "%X"))

##get PRS



##compute predictions for all data
print("predict PGS for all")
print(format(Sys.time(), "%X"))
pred <- big_prodVec(G_imp, info_snps$beta_auto,  ind.col = info_snps[["_NUM_ID_"]],ncores = opt$threads)

pred2 <- as.data.frame(pred) %>%
                mutate(cpr_enc=obj.bigSNP$fam[,"sample.ID"])




print("Done")
print(format(Sys.time(), "%X"))
saveRDS(pred2,opt$out_pred)

