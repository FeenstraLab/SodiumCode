## format GWAS summary statistics
## Matching variants between genotype data and summary statistics
## select only hm3 variants

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Starting sumstats parsing!\n"))

library(optparse)
library(stringr)
library(dplyr)
library(tidyr)
library(bigsnpr)
library(bigreadr)
library(ggplot2)
library(gridExtra)

parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character',
              help = "GWAS summary statistics file path"),
  make_option("--bed",
              action = "store",
              default = NA,
              type = 'character',
              help = "genotype data with hm3 variants directory"),
  make_option("--freq",
              action = "store",
              default = NA,
              type = 'character',
              help = "minor allele frequency"),
  make_option("--rds_out",
              action = "store",
              default = NA,
              type = 'character',
              help = "output bigsnor format genotype "),
    make_option(c("-t", "--threads"),
              action = "store",
              default = 1,
              type = 'integer',
              help = "Number of threads to use"),
    make_option("--QC_plot",
              action = "store",
              default = NA,
              type = 'character',
              help = "output QC plot"),
  make_option("--out_noQC",
              action = "store",
              default = NA,
              type = 'character',
              help = "output sumstats"),
  make_option("--out_QC1",
              action = "store",
              default = NA,
              type = 'character',
              help = "output sumstats"),
  make_option("--out_QC2",
              action = "store",
              default = NA,
              type = 'character',
              help = "output sumstats")
)

if (interactive()) {

} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Reading genotype file, impute, then calucated allle frequency, output map file with minor allele frequency  ..."))


## input genotype data,save to RDS format
rds <- snp_readBed2(opt$bed,
               backingfile = opt$rds_out,
               ncores = opt$threads)

obj.bigSNP <- snp_attach(rds)
map <- obj.bigSNP$map %>%
        transmute(chr = as.integer(chromosome), pos = physical.pos, a1 = allele1, a0 = allele2,SNP=marker.ID)

# # imputation
# G <- obj.bigSNP$genotypes
# G2 <- snp_fastImputeSimple(G,ncores = opt$threads)
# map$freq <- snp_MAF(G2,ncores=opt$threads)
freq <- read.table(opt$freq,header=T) %>%
        select(SNP,MAF)

MAP <- merge(map,freq,by="SNP")
# #for (i in 1:2){
# for (i in 1:22){
#     obj.bigSNP <- snp_attach(paste0(opt$genotype_hm3,"/chr",i,".rds"))

#     map <- obj.bigSNP$map %>%
#         transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,)

#     # imputation
#     G <- obj.bigSNP$genotypes
#     G2 <- snp_fastImputeSimple(G,ncores = opt$threads)

#     # maf
#     map$freq <- snp_MAF(G2,ncores=opt$threads)

#     if(i==1){
#       MAP <- map
#     }else{
#       MAP <- rbind(MAP,map)
#     }

# }



## input GWAS summary data
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Reading sumstats input file ..."))

#gwas <- read.table("tmp/gwas.txt.gz",header=T) %>%
gwas <- fread2(opt$gwas,
      select=c("MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","TotalSampleSize","HetPVal"),
      col.names=c("MarkerName","a1","a0","freq","beta","beta_se","N","HetPVal")
        ) %>%
        separate(MarkerName,c("CHR","BP")) %>%
        mutate(across('CHR', str_replace, 'chr', '')) %>%
        mutate(across('a0',toupper)) %>%
        mutate(across('a1',toupper)) %>%
        mutate(chr=as.numeric(CHR),pos=as.numeric(BP)) %>%
        filter(chr %in% 1:22,pmin(freq,1-freq) >0.01) %>%
        select(chr,pos,a1,a0,beta,beta_se,N,HetPVal) 

# ## To get rsid for sumstats
# hm3 <- readRDS(opt$hm3) %>%
#   rename(pos_hg19=pos) %>%
#   rename(pos=pos_hg38) %>%
#   select(chr,pos,a0,a1,rsid)

# gwas2 <- snp_match(gwas,hm3)

cat(" Complete!\nPreview of input file:\n")
head(gwas)




## Matching variants between genotype data and summary statistics
## https://github.com/privefl/paper-misspec/blob/main/code/prepare-sumstats/vitaminD.R

info_snp <- as_tibble(bigsnpr::snp_match(gwas, MAP))

png("results/d20230804/PGS/parsed_sumstats/excl3_lm_le50_weights/hist_N.png")
hist(info_snp$N, breaks = 50, main = NULL, xlab = "Sample size")
abline(v = 0.7 * max(info_snp$N), col = "red")
abline(v = 0.75 * max(info_snp$N), col = "blue")
dev.off()


info_snp2 <- info_snp %>%
  filter(HetPVal>1e-5) %>%
  mutate(n_eff=N,
         sd_af = sqrt(2 * MAF * (1 - MAF)),
         sd_ss = 1 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5))

info_snp2$is_bad <- with(info_snp2,
                         sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                        sd_ss < 0.1 | sd_af < 0.05)

# info_snp2$is_bad2 <- with(info_snp2,
#                          sd_ss < (0.5 * sd_af) | sd_ss > (sd_af + 0.1) |
#                         sd_ss < 0.1 | sd_af < 0.05)

p1 <- qplot(sd_af, sd_ss,  color = ifelse(is_bad, "Yes", "No"),alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr(0.8) +
  coord_equal() +
  theme(legend.position = c(0.2, 0.8)) +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

# p2 <- qplot(sd_af, sd_ss,  color = ifelse(is_bad2, "Yes", "No"),alpha = I(0.5),
#       data = slice_sample(info_snp2, n = 100e3)) +
#   theme_bigstatsr(0.8) +
#   coord_equal() +
#   theme(legend.position = c(0.2, 0.8)) +
#   scale_color_viridis_d(direction = -1) +
#   geom_abline(linetype = 2, color = "red") +
#   labs(x = "Standard deviations derived from allele frequencies",
#        y = "Standard deviations derived from the summary statistics",
#        color = "Removed?")

# p3 <- grid.arrange(p1,p2, nrow = 1,ncol=2)
 
ggsave(opt$QC_plot, p1, width = 7, height = 7)

saveRDS(info_snp2, opt$out_noQC)
saveRDS(filter(info_snp2, !is_bad), opt$out_QC1)
saveRDS(filter(info_snp2, !is_bad, N > (0.7 * max(N))), opt$out_QC2)

