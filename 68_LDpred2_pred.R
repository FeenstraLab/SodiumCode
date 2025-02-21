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
  make_option("--sample",
              action = "store",
              default = NA,
              type = 'character',
              help = 'sample used for estimate LD relationship matrix'),
  make_option(c("-t", "--threads"),
              action = "store",
              default = 1,
              type = 'integer',
              help = "Number of threads to use"),
  make_option("--phen",
              action = "store",
              default = 1,
              type = 'character',
              help = "phenotype data with sodium info"),
  make_option("--sumstats",
              action = "store",
              default = 1,
              type = 'character',
              help = "Number of threads to use"),
  make_option("--auto",
              action = "store",
              default = NA,
              type = 'character',
              help = 'input of LDpred2 auto model'),
  make_option("--out_beta",
              action = "store",
              default = NA,
              type = 'character',
              help = 'output of LDpred2 auto beta effect'),
  make_option("--out_pred",
              action = "store",
              default = NA,
              type = 'character',
              help = ' PRS output'),
  make_option("--quantile_plot",
              action = "store",
              default = NA,
              type = 'character',
              help = 'quantile plot for validation dataset in samples ')
  
)

if (interactive()) {

} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

## input matrix
obj.bigSNP <- snp_attach(opt$genotype_RDS)
G <- obj.bigSNP$genotypes

## input auto mdoel results
auto <-  readRDS(opt$auto)
beta_auto <- sapply(auto, function(auto) auto$beta_est)


## define validation and test data set
## sample with sodium records can be validation 
## sample missing sodium records can be test
phen <- fread2(opt$phen)
sample <- fread2(opt$sample)

sodium_mean <- phen %>%
          filter(cpr_enc %in% sample$cpr_enc) %>%
          select(cpr_enc,VALUE) %>%
          group_by(cpr_enc) %>%
          summarise(sodium=mean(VALUE))
sodium_cov <- phen %>%
          filter(cpr_enc %in% sample$cpr_enc) %>%
          select(cpr_enc,sex,sodium_age,paste0("PC",1:10)) %>%
          group_by(cpr_enc) %>%
          arrange(sodium_age) %>%
          filter(row_number()==1) %>%
          inner_join(y=sodium_mean,by="cpr_enc") %>%
          select(!sex)

# sample <- sample %>%
#           left_join(y=sodium_cov,by="cpr_enc") 


fam <- obj.bigSNP$fam %>%
        mutate(order=1:nrow(obj.bigSNP$fam)) %>%
        left_join(y=sodium_cov,join_by(sample.ID==cpr_enc)) %>%
        arrange(order)

ind.val <- fam[!is.na(fam$sodium),"order"]
ind.test <- fam[is.na(fam$sodium),"order"]


##
sumstats <- readRDS(opt$sumstats) %>%
      rename(rsid=SNP) %>%
      select("chr","pos", "rsid","a0", "a1","beta", "beta_se","n_eff")
map <- obj.bigSNP$map %>%
        rename(chr=chromosome,rsid=marker.ID,pos=physical.pos,a1=allele1,a0=allele2)

info_snps <- snp_match(sumstats, map,join_by_pos = FALSE)


## impute G
print("impute G")
print(format(Sys.time(), "%X"))
G_imp <- snp_fastImputeSimple(
  G,
  ncores = opt$threads
)
print(format(Sys.time(), "%X"))

##get PRS
print("predict PGS for validation")
print(format(Sys.time(), "%X"))
pred_auto <- big_prodMat(G_imp, beta_auto,ind.row = ind.val,ind.col = info_snps[["_NUM_ID_"]],ncores = opt$threads)
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- as.data.frame(rowMeans(beta_auto[, keep])) %>%
        bind_cols(info_snps) %>%
        select(!c("pos.ss","_NUM_ID_.ss","_NUM_ID_"))
names(final_beta_auto)[1] <- "beta_auto"
saveRDS(final_beta_auto,opt$out_beta)




print("accuray of PGS")
print(format(Sys.time(), "%X"))
pred_auto_val <- big_prodVec(G_imp, final_beta_auto$beta_auto, ind.row = ind.val, ind.col = info_snps[["_NUM_ID_"]],ncores = opt$threads)
print("Correlation between PRS and sodium mean value")
print(pcor(pred_auto_val,fam$sodium[ind.val],NULL))
#0.1095197 0.0794343 0.1394058
print(pcor(pred_auto_val,fam$sodium[ind.val],fam[ind.val,c("sex","sodium_age",paste0("PC",1:10))]))
#[1] 0.11096941 0.08084888 0.14088745




##compute predictions for all data
print("predict PGS for all")
print(format(Sys.time(), "%X"))
pred_auto_all <- big_prodVec(G_imp, final_beta_auto$beta_auto,  ind.col = info_snps[["_NUM_ID_"]],ncores = opt$threads)

pred_auto_all2 <- as.data.frame(pred_auto_all) %>%
                mutate(cpr_enc=obj.bigSNP$fam[,"sample.ID"]) %>%
                mutate(sex=obj.bigSNP$fam[,"sex"]) 



print("Done")
print(format(Sys.time(), "%X"))
write.table(pred_auto_all2,opt$out_pred,col.names=T,row.names=F,quote=F,sep="\t")


png("tmp/beta.png")
hist(pred_auto_all2$pred_auto_all)
dev.off()

## quantile plot
get_quantile <- function(x, num.quant, quant.ref){
  quant <- as.numeric(cut(x,
                          breaks = unique(quantile(
                            x, probs = seq(0, 1, 1 / num.quant)
                          )),
                          include.lowest = T))
  if(is.null(quant.ref)){
    quant.ref <- ceiling(num.quant / 2)       #define reference factor
  }
  quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
  return(quant)
}

dat <- merge(pred_auto_all2,sodium_cov,by="cpr_enc")
dat$quantile <- get_quantile(dat[,"pred_auto_all"],num.quant=10, quant.ref=1)
dat$age2 <- dat$sodium_age^2
lr <- lm(sodium ~ sex +sodium_age+age2+PC1+PC2+PC3+PC4+PC5,data=dat)
dat$lr_zscore <- (resid(lr)-mean(resid(lr)))/sd(resid(lr))

lr <- summary(lm(lr_zscore ~ -1 + quantile ,data=dat))
lr_effect <- lr$coefficients[,1]
lr_sd <- lr$coefficients[,2]
q975 <- qnorm(0.975)
s <- cbind(c(1:10),lr_effect,lr_sd,lr_effect+q975*lr_sd,lr_effect-q975*lr_sd)
s <- as.data.frame(s)
names(s) <- c("Quantile","Effect","SD","CI.U","CI.L")

p1 <- ggplot(s,aes(x=Quantile,y=Effect,ymin=CI.L,ymax=CI.U))+
  geom_point(position=position_dodge(0.5))+
  geom_errorbar(width=0.1,position=position_dodge(0.5)) +
  xlab("PGS Quantiles")+
  ylab("Effect (95% CI)")+
  scale_x_continuous(breaks=seq(1,10,1)) +
  theme(legend.title = element_blank())

ggsave(opt$quantile_plot,p1,width=4,height=4)