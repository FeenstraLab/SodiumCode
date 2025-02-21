print("Remove multiallelic loci and indel loci which has either only 1 line or more than 2 lines!")
print("save SNPs which has 2 lines and merge them!")
print("Manhattan plot")

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Initiating resolution of multiallelic sites...\n"))

### Packages and options -------------------------------------------------------
library(optparse)
library(bigreadr)
library(dplyr)
library(tidyr)
library(R.utils)
library(qqman)


parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--gwas_padj",
              action = "store",
              default = NA,
              type = 'character'),
    make_option("--manplot",
              action = "store",
              default = NA,
              type = 'character',
              help = "manhattan plot"),
    make_option("--qqplot",
              action = "store",
              default = NA,
              type = 'character',
              help = "qq plot ")
)

if (interactive()) {
  
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

### Import and Loading ---------------------------------------------------------
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Reading input file ", opt$gwas," ..."))

df <-  fread2(opt$gwas)

cat("Complete!\nPreview of input file:\n")

#head(df)

cat(paste0("Your initial dataset has ", nrow(df), " variants.\n"))

### Resolution of multiallelic sites -------------------------------------------
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Tagging variants...\n"))

# First, tag each variant 
# :IG: multiallelic insertion   # for n=2, merge; if n!=2, remove it
# :IG.0: or :IG.1: or :IG.2: multallelic deletion     # for n=2, merge; if n!=2, remove it
# :SG: multiallelic             # for n=2, merge; if n!=2, remove it
# :SG      biallelic                Keep all
# :IG or :IG.0   indel(insertion/deletion), Keep all

df <-  df %>%
  filter(A1FREQ>=0.01 &  A1FREQ<=0.99) %>%
  mutate(MARKER = case_when(grepl(":SG:|:IG:|:IG.0:|:IG.1:|:IG.2:", ID) ~ "multiallelic",     #including multiallelic also for insertion and deletion
                            TRUE ~ "biallelic"))    ## including biallelic and indel with only 1 record
#  group_by(CHROM, GENPOS) %>%
#  mutate(n=n())
cat(paste0("After remove maf<0.01 dataset has ", nrow(df), " variants.\n"))


(table(df$MARKER))
# biallelic multiallelic 
#    10999658      6098571  

#table(df[df$MARKER=="biallelic","n"])
#table(df[df$MARKER=="indel","n"])

# Within the "multiallelic" category we have: 
# 1. Sites for which multiple alleles are tested in the GWAS.
# 2. Sites for which very low frequency alleles have been filtered out & only 2 
#    are tested in the GWAS. There are two alternatives here:
#    a. Either the removed alleles had a negligible effect & the two remaining 
#       capture most of the effect.
#    b. Or the removed alleles had a considerable effect. This can be observed in
#       the summary statistics because the effect sizes of the leftover alleles
#       are not mirror images of each other.

# (2a) is, in practical terms, biallelic and can be resolved as explained in the 
# README

# Second, resolve as far as possible
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Resolving multiallelic...\n"))

resolved_multiallelic = df %>%
  filter(MARKER == "multiallelic") %>% # get multiallelic sites
#  slice_head(n=100000) %>%
  group_by(CHROM, GENPOS) %>%
  mutate(n=n()) %>% # how many alleles do we observe?
  filter(n==2) %>% # get only observed biallelic 
  mutate(ALLELE0 = case_when(!is.na(lead(ALLELE1)) ~ lead(ALLELE1),
                             !is.na(lag(ALLELE1)) ~ lag(ALLELE1))) %>% # recode
  # we are setting the max difference in BETA to 1.2 fold. The rest, we don't resolve. 
  mutate(diff_BETA = abs(BETA/lead(BETA)),
         diff_BETA = case_when(is.na(diff_BETA) ~ lag(diff_BETA),
                               TRUE ~ diff_BETA)) %>%
  filter((INFO >= 0.8) & (diff_BETA <= 1.2) & (diff_BETA >= (1/1.2))) %>% # remove what we can't resolve
  separate(ID,c("CHR","BP","site","site_n"),sep=":") %>%
  mutate(ID = paste("chr", CHROM, ":", GENPOS, ":",site, sep = "")) %>% # rename the SNP so it looks like a biallelic site
  distinct(CHROM, GENPOS, .keep_all = TRUE) %>% # remove duplicates
  select(-c(n, diff_BETA,CHR,BP,site,site_n,EXTRA))

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Succesfully resolved a total of multiallelic ", nrow(resolved_multiallelic), " sites.\n"))


# Merge resolved sites back in
# BF 240223: Also removing instances where one allele is "*"
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Creating final dataset...\n"))

clean_df = df %>% 
  filter(MARKER=="biallelic") %>%
  select(-c(EXTRA)) %>%
 # anti_join(resolved_multiallelic, by = c("CHROM", "GENPOS", "MARKER")) %>% # what didn't need resolving/could not be resolved
  rbind(resolved_multiallelic) %>% # merge resolved multiallelic sites back in
  filter(!is.na(ALLELE0) & !is.na(ALLELE1)) %>%  #18  chr18:80256321:IG 0 80256321  N NA
  filter(ALLELE0!="*" & ALLELE1!="*") %>%  # BF added on 240223
  select(-MARKER) %>% # remove unwanted variables
  arrange(CHROM, GENPOS) # sort


print("Check if any ! left")
anymark <- clean_df[clean_df$ALLELE0=="!" | clean_df$ALLELE1=="!",]
(nrow(anymark))
(head(anymark))

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Your final dataset has ", nrow(clean_df), " variants.\n"))



### 
gwas <- clean_df %>%
 #       mutate(SNP=ID) %>%
        rename(SNP=ID,CHR=CHROM,BP=GENPOS) %>%
        select(SNP,CHR,BP,BETA,SE,P,N,ALLELE0,ALLELE1,A1FREQ,INFO)
 
##get adjusted p value and se value
adj <- function(data,pcol,betacol,secol){
        p_quant <- qchisq(data[,pcol],1,low=F)
        lambda <- median(p_quant)/qchisq(0.5,1)
  if(lambda>1){
      data$p_adj <- pchisq(p_quant/lambda,1,low=F)
      data$se_adj <- abs(data[,betacol]/qnorm(data$p_adj/2,low=F))
  }else{
      data$p_adj <- data[,pcol]
      data$se_adj <- data[,secol]
  }
        return(list(lambda,data))
}

gwas_padj <- adj(data=gwas,pcol="P",betacol="BETA",secol="SE")[[2]]

lambda <- adj(data=gwas,pcol="P",betacol="BETA",secol="SE")[[1]]
print(lambda)



### Writing output -------------------------------------------------------------
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Writing output ", opt$gwas," ..."))

write.table(gwas_padj, opt$gwas_padj,col.names=T,row.names=F,quote=F,sep="\t")
gzip(opt$gwas_padj)

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "All done!"))



## manhattan plot
print("plot for P")
maxy <- max(-log10(gwas$P),na.rm=T)+1
png(file=opt$manplot, width =720, height = 480)
manhattan(gwas,ylim=c(0,maxy),col= c("blue4", "orange3"))
dev.off()

png(file=opt$qqplot, width = 480, height = 480)
qq(gwas$P, main = "Q-Q plot of GWAS p-values",col = "blue4")
dev.off()