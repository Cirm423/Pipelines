log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(EnhancedVolcano)
library(glue)

res <- read.table(snakemake@input[[1]], sep='\t', header=1, )

title = glue("{snakemake@params[['contrast']][1]}-vs-{snakemake@params[['contrast']][2]}")

svg(snakemake@output[[1]])
EnhancedVolcano(res, lab = res$gene, title = title, x = 'log2FoldChange', y= 'padj',pCutoff = 0.05,FCcutoff=1.0)
dev.off()
