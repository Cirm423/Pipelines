log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)

res <- read.table(snakemake@input[["res"]], sep='\t', header=1, )

#Remove version from ensembl IDs
res$gene <- sub('\\.[0-9]*$', '', res$gene)

gtf <- readr::read_tsv(snakemake@input[["gtf"]], col_names=FALSE, comment="#", show_col_types=FALSE)

# Keep only gene features
genes <- gtf[gtf[[3]] == "gene", ]

# Extract attribute values from the 9th column
extract_attr <- function(x, key) {
  pat <- paste0(key, " \"([^\"]+)\"")
  m <- regmatches(x, regexec(pat, x))
  vapply(m, function(y) if (length(y) >= 2) y[2] else NA_character_, character(1))
}

gene_map <- data.frame(
  ensembl_gene_id = sub("\\.[0-9]+$", "", extract_attr(genes[[9]], "gene_id")),
  hgnc_symbol     = extract_attr(genes[[9]], "gene_name"),
  stringsAsFactors = FALSE
)

# Keep one symbol per Ensembl ID and drop empty symbols
gene_map <- gene_map %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

res_complete <- left_join(res, gene_map, by = join_by(gene == ensembl_gene_id))

# Convert empty strings to NA
res_complete <- res_complete %>% mutate_if(is.character, list(~na_if(.,"")))

# Fall back to Ensembl ID when symbol is missing
res_complete$hgnc_symbol[is.na(res_complete$hgnc_symbol)] <- res_complete$gene[is.na(res_complete$hgnc_symbol)]

res$gene <- res_complete$hgnc_symbol

write.table(res, file=snakemake@output[[1]], row.names=FALSE, sep='\t')