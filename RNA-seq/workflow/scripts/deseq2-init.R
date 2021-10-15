log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("tximport")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)
files <- snakemake@input[["counts"]]
names(files) <- basename(dirname(files))
#Making sure both coldata and files are ordered in the same way
files <- files[order(match(names(files),row.names(coldata)))]
#cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
#Converst genes with length 0 to length 1 to make the matrix, they will be removed downstream anyway. Makes sure dds creation does not fail.
txi.rsem$length[txi.rsem$length == 0] <- 1
# dds <- DESeqDataSetFromMatrix(countData=cts,
#                               colData=coldata,
#                               design=as.formula(snakemake@params[["model"]]))

dds <- DESeqDataSetFromTximport(txi.rsem,
                                colData=coldata,
                                design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

# Write dds object as RDS
saveRDS(dds, file=snakemake@output[[1]])
# Write normalized counts
norm_counts = counts(dds, normalized=T)
write.table(data.frame("gene"=rownames(norm_counts), norm_counts), file=snakemake@output[[2]], sep='\t', row.names=F)
