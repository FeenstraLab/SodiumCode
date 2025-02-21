cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Initiating liftover...\n"))

### Packages and options ---------------
library(optparse)
library(rtracklayer)
library(dplyr)
library(stringr)

parsed_opts <- list(
  make_option("--input",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--chain_path",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--chr_lookup",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--output",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--discarded",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--discard_info",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--from_build",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--to_build",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--chr_column",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--pos_column",
              action = "store",
              default = NA,
              type = 'character')
)

if (interactive()) {

} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Debugging info:\nopt <- "))
dput(opt)

if (opt$from_build == "hg37") opt$from_build <- "hg19"
if (opt$to_build == "hg37") opt$to_build <- "hg19"


### Import and Loading ---------------
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Reading input file ", opt$input," ..."))

raw_input <- data.table::fread(opt$input, data.table = F)
raw_input$CHR <- apply(as.data.frame(raw_input$MarkerName),1,function(x){unlist(strsplit(x,":"))[1]})
raw_input$BP <- apply(as.data.frame(raw_input$MarkerName),1,function(x){unlist(strsplit(x,":"))[2]})

cat("Complete!\nPreview of input file:\n")

head(raw_input)

raw_input_columns <- colnames(raw_input)
raw_without_pos <- raw_input_columns[-which(raw_input_columns == opt$pos_col)]
unique_out_columns <- make.names(c("seqnames", "start", "end", "width", "strand", raw_without_pos),
                                 unique = T)

processed_column_order <- c(opt$chr_column, 
                            opt$pos_column, 
                            unique_out_columns[6:length(unique_out_columns)])

raw_input$input_row_index <- seq_len(nrow(raw_input))

### Parse chromosome names
chr_lookup <- readr::read_tsv(opt$chr_lookup)%>% 
  rbind(data.frame(input_chr = "chr23", internal_liftover_parsed_chr = "chrX")) ## NIAN addition 240306 to handle chrX
colnames(chr_lookup)[1] <- opt$chr_column

raw_input[, opt$chr_column] <- as.character(raw_input[, opt$chr_column])

parsed_chr_input <- raw_input %>% 
  left_join(chr_lookup)

ready_for_liftover <- parsed_chr_input %>% 
  filter(! is.na(internal_liftover_parsed_chr))

first_discarded_rows <- parsed_chr_input %>% 
  filter(is.na(internal_liftover_parsed_chr)) %>% 
  mutate(discard_reason = "invalid chromosome")

if (nrow(first_discarded_rows) > 0) {
  cat(paste0("[", format(Sys.time(), "%X"), "] ",
             "NOTE: Discarded ", nrow(first_discarded_rows), 
             " variants during chromosome parsing (", 
             round(nrow(first_discarded_rows) / nrow(raw_input) * 100, 2) ,
             "% of the ", nrow(raw_input), " original variants in the input file).\n"))
}

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Running liftover..."))

liftover_chain <- import.chain(opt$chain_path)

input_granges <- 
  GenomicRanges::makeGRangesFromDataFrame(ready_for_liftover,
                                          seqnames.field = "internal_liftover_parsed_chr",
                                          start.field = opt$pos_column,
                                          end.field = opt$pos_column,
                                          ignore.strand = TRUE,
                                          keep.extra.columns = TRUE)

gr_lifted <- liftOver(x = input_granges, chain = liftover_chain) %>% 
  unlist() %>% 
  data.frame()

discarded_index <- setdiff(ready_for_liftover$input_row_index, gr_lifted$input_row_index)

second_discarded_rows <- ready_for_liftover %>% 
  filter(input_row_index %in% discarded_index) %>%
  mutate(discard_reason = "removed by liftover function")

all_discards <- bind_rows(first_discarded_rows, second_discarded_rows)

cat("Complete!\n")

if (length(discarded_index) > 0) {
  cat(paste0("\n[", format(Sys.time(), "%X"), "] ",
             "NOTE: Discarded ", length(discarded_index), 
             " variants during liftover step (", 
             round(length(discarded_index) / nrow(raw_input) * 100, 2) ,
             "% of the ", nrow(raw_input), " original variants in the input file).\n"))
}

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Writing output files..."))

colnames(gr_lifted)[1] <- opt$chr_column
colnames(gr_lifted)[2] <- opt$pos_column

gr_lifted %>% 
  select(all_of(processed_column_order)) %>% 
  readr::write_tsv(file = opt$output)

if (nrow(all_discards) > 0) {
  readr::write_tsv(all_discards, file = opt$discarded)
  
  cat(paste0("[", Sys.time(), "] ",
             "Info - discarded variants during liftover.",
             "\n", paste0(rep("-", 80), collapse = ""),
             "\ninput file: ", opt$input, 
             "\nfrom: ", opt$from_build,
             "\nto: ", opt$to_build,
             "\noutput file: ", opt$output, 
             "\ndiscarded rows in file: ", opt$discarded,
             "\n\ntotal variants discarded: ", nrow(all_discards), 
             " (", round(nrow(all_discards) / nrow(raw_input) * 100, 2), 
             "% of ", nrow(raw_input), " variants in the input file).", 
             "\n", paste0(rep("-", 80), collapse = "")), 
      file = opt$discard_info)
  
  cat("\nBreakdown of discarded variants by input chromosome:\n\n", 
      file = opt$discard_info, append = T)
  readr::write_tsv(count(all_discards, .data[[opt$chr_column]], sort = T, name = "count"), 
                   file = opt$discard_info, append = T, col_names = T)
  
  cat(paste0(rep("-", 80), collapse = ""), 
      file = opt$discard_info, append = T)
  cat("\nBreakdown of discarded variants by input chromosome AND discard reason:\n\n", 
      file = opt$discard_info, append = T)
  readr::write_tsv(count(all_discards, .data[[opt$chr_column]], discard_reason, sort = T, name = "count"), 
                   file = opt$discard_info, append = T, col_names = T)
} else {
  file.create(opt$discarded)
  file.create(opt$discard_info)
}

cat("Complete!\n")
