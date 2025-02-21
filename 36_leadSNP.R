##==============================leadSNP on each chromosome
library(bigreadr)
library(dplyr)
library(tidyr)
library(optparse)
parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--leadSNP",
              action = "store",
              default = NA,
              type = 'character')
)

if (interactive()) {

} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}




gwas <- fread2(opt$gwas)  %>%
  rename(P="P-value") %>%
  filter(P < 5e-8) %>%
  separate(MarkerName,c("CHR","BP"),sep=":") %>%
  mutate(BP=as.numeric(BP)) %>%
  arrange(CHR,BP) %>%
  group_by(CHR) %>%
  mutate(distance=(BP - lag(BP)))

# Genome-wide significant loci were defined as regions with one or more SNPs with p-value <5×10−8, and these SNPs were defined as belonging to different loci if the distance between them was >500 kb
indexfun <- function(sdat){
  sdat <- sdat[order(sdat$BP),]
  g <- 1
  for (r in 1:nrow(sdat)){
    if(is.na(sdat[r,"distance"])){
      sdat[r,"index"] <- g
    }else if(sdat[r,"distance"]<=500000){
      sdat[r,"index"] <- g
    }else if(sdat[r,"distance"]>500000){
      g <- g+1
      sdat[r,"index"] <- g
    }
  }
   return(sdat)
} 
 
df2 <- gwas %>%
    group_by(CHR) %>%
    do(indexfun(sdat=.))

top <- df2 %>%
  group_by(CHR,index) %>%
  slice(which.min(P))


#print(as.data.frame(top))
write.table(top,opt$leadSNP,col.names=T,row.names=F,quote=F)
