
## Calculate LD Matrix
library(optparse)
library(bigsnpr)
library(dplyr)
library(readr)
library(stringr)

parsed_opts <- list(
  make_option("--genotype_RDS",
              action = "store",
              default = NA,
              type = 'character',
              help = 'RDS file path for genotype data of a single chromosome'),
  make_option("--matrix_out",
              action = "store",
              default = NA,
              type = 'character',
              help = 'RDS file path for storing LD matrix'),
  make_option("--genetic_pos",
              action = "store",
              default = NA,
              type = 'character',
              help = 'Directory containing the Snakefile (reference for other project files).'),
  make_option(c("-t", "--threads"),
              action = "store",
              default = 1,
              type = 'integer',
              help = "Number of threads to use"),
   make_option("--chr",
              action = "store",
              default = 1,
              type = 'integer',
              help = "chromosome"),
  make_option("--sumstats",
              action = "store",
              default = 1,
              type = 'integer',
              help = "Number of threads to use")
)

if (interactive()) {

} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

### Data preparation ------------
obj.bigSNP <- snp_attach(opt$genotype_RDS)

#' Isn't this fascinating?
#+ fasc2, echo=TRUE
G   <- obj.bigSNP$genotypes
CHR <- as.integer(obj.bigSNP$map$chromosome)
POS <- obj.bigSNP$map$physical.pos
rsid <- obj.bigSNP$map$marker.ID

map <- obj.bigSNP$map %>%
      rename(chr=chromosome,rsid=marker.ID,pos=physical.pos,a0=allele2,a1=allele1)


## sumstats
sumstats <- readRDS(opt$sumstats) %>%
      rename(rsid=SNP) %>%
      select("chr","pos", "rsid","a0", "a1","beta", "beta_se")
info_snps <- snp_match(sumstats, map,join_by_pos = FALSE)
#_NUM_ID_.ss index in sumstats
#_NUM_ID_ index in map

#matching 
POS2 <- snp_asGeneticPos(CHR,POS, dir =opt$genetic_pos)

### Calculate correlation matrix --------

#for (chr in 1:22){
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Calculating LD matrix..."))

  ind.chr <- which(info_snps$chr==opt$chr)
  ind.chr2 <- info_snps$`_NUM_ID_`[ind.chr]
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = opt$threads)

#   if (chr == 1) {
#       ld <- Matrix::colSums(corr0^2)
#       corr <- as_SFBM(corr0, tmp, compact = TRUE)
#     } else {
#       ld <- c(ld, Matrix::colSums(corr0^2))
#       corr$add_columns(corr0, nrow(corr))
#     }
# #}

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Saving LD matrix to ", opt$matrix_out, "..."))

saveRDS(corr0, file = opt$matrix_out)

