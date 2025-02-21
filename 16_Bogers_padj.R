print("Adjust p value for Bogers' GWAS")

library(optparse)
library(qqman)
library(R.utils)

parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character',
              help = "Regenie GWAS summary file"),
   make_option("--gwas_padj",
              action = "store",
              default = NA,
              type = 'character',
              help = "output with adjusted p value and se value"),
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
    opt <- list(gwas="results/d20230523/Bogers/Osmolality_lifted_to_hg38.tsv.gz",
    		gwas_padj="results/d20230523/Bogers/gwas_padj.txt",
        manplot="results/d20230523/Bogers/manplot.png",
        qqplot="results/d20230523/Bogers/qqplot.png"
    )
  } else {
    opt <- parse_args(OptionParser(option_list = parsed_opts))
    
}

#read file
gwas <- read.table(gzfile(opt$gwas),header=T)
names(gwas)[2] <- "POS_hg38"



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


write.table(gwas_padj,opt$gwas_padj,col.names=T,row.names=F,quote=F)
gzip(opt$gwas_padj)


################
print("plot for P")

dat <- gwas[,c("SNP","CHR","POS_hg38","P")]
names(dat) <- c("SNP","chr","BP","P")

## change "chr2" to "2"
chr <- as.data.frame(matrix(NA,22,2))
names(chr) <- c("chr","CHR")
chr$chr <- paste0("chr",1:22)
chr$CHR <- 1:22
dat <- merge(dat,chr,by="chr")

maxy <- max(-log10(dat$P),na.rm=T)+1
png(file=opt$manplot, width =720, height = 480)
manhattan(dat,ylim=c(0,maxy),col= c("blue4", "orange3"))
dev.off()

png(file=opt$qqplot, width = 480, height = 480)
qq(dat$P, main = "Q-Q plot of GWAS p-values",col = "blue4")
dev.off()


