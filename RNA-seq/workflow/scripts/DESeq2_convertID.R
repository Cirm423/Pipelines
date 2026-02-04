log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(biomaRt)
library(dplyr)

res <- read.table(snakemake@input[[1]], sep='\t', header=1, )

#Remove version from ensembl IDs
res$gene <- sub('\\.[0-9]*$', '', res$gene)

human_mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"  # Human dataset
)

gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res$gene,
  mart = human_mart,
  uniqueRows=TRUE
)

gene_symbols_unique <- gene_symbols %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

res_complete <- left_join(res, gene_symbols_unique, by=join_by(gene == ensembl_gene_id))

#Remove IDs that are empty strings
res_complete <- res_complete %>% mutate_if(is.character, list(~na_if(.,""))) 

#Paste EnsemblIDs on missing symbols
res_complete$hgnc_symbol[is.na(res_complete$hgnc_symbol)] <- res_complete$gene[is.na(res_complete$hgnc_symbol)]

res$gene <- res_complete$hgnc_symbol

write.table(res, file=snakemake@output[[1]], row.names=FALSE, sep='\t')