library(optparse)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bigreadr)

parsed_opts <- list(
  make_option("--genotype_RDS",
              action = "store",
              default = NA,
              type = 'character',
              help = 'RDS file path for genotype data of a single chromosome'),
  make_option("--matrix_dir",
              action = "store",
              default = NA,
              type = 'character',
              help = 'RDS file path for storing LD matrix'),
  make_option(c("-t", "--threads"),
              action = "store",
              default = 1,
              type = 'integer',
              help = "Number of threads to use"),
  make_option("--sumstats",
              action = "store",
              default = 1,
              type = 'character',
              help = "filtered summary data"),
  make_option("--out_auto",
              action = "store",
              default = NA,
              type = 'character',
              help = 'output of LDpred2 auto model'),
    make_option("--plot_chain",
              action = "store",
              default = NA,
              type = 'character',
              help = ' verify whether the chains “converged” by looking at the path of the chains')
  
)

if (interactive()) {


} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}


## merge ld matrix
for(chr in 1:22){
  corr0 <- readRDS(paste0(opt$matrix_dir,"/chr",chr,".rds"))

  if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, opt$matrix_dir, compact = TRUE)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }

}

print(format(Sys.time(), "%X"))

sumstats <- readRDS(opt$sumstats) %>%
      rename(rsid=SNP) %>%
      select("chr","pos", "rsid","a0", "a1","beta", "beta_se","n_eff")

obj.bigSNP <- snp_attach(opt$genotype_RDS)
map <- obj.bigSNP$map %>%
        rename(chr=chromosome,rsid=marker.ID,pos=physical.pos,a1=allele1,a0=allele2)

info_snps <- snp_match(sumstats, map,join_by_pos = FALSE)

## Estimate of h2 from LD Score regression
(ldsc <- with(info_snps, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))

ldsc_h2_est <- ldsc[["h2"]]
print(format(Sys.time(), "%X"))

#===================automatic model, it takes more than 1 hour with burn in 500, iter 200
coef_shrink <- 0.95
set.seed(1)  
multi_auto <- snp_ldpred2_auto(
  corr, info_snps, h2_init = ldsc_h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = 30), ncores = opt$threads,
  # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  burn_in = 800, num_iter = 500, report_step = 20,
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)
#str(multi_auto, max.level = 1)

print(format(Sys.time(), "%X"))

saveRDS(multi_auto,opt$out_auto)


##================================== verify results

auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
ggsave(opt$plot_chain,width=7,height=4)


#an automatic way of filtering bad chains by comparing the scale of the resulting predictions. We have tested a somewhat equivalent and simpler alternative since, which we recommend here:
#(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
#(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))



#To get the final effects / predictions, you should only use chains that pass this filtering:
#beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))


#saveRDS(beta_auto,opt$out_betas)

# #Parameters h2, p, and α  (and 95% CIs) can be estimated using
# all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
# quantile(all_h2, c(0.5, 0.025, 0.975))
# all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
# quantile(all_p, c(0.5, 0.025, 0.975))

# all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
# quantile(all_alpha, c(0.5, 0.025, 0.975))

# Predictive performance r2
#  can also be inferred from the Gibbs sampler:

# bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
# all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
#   b1 <- bsamp[[ic]]
#   Rb1 <- apply(b1, 2, function(x)
#     coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
#   b2 <- do.call("cbind", bsamp[-ic])
#   b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
# }))
# quantile(all_r2, c(0.5, 0.025, 0.975))

