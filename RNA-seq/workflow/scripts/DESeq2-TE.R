#! /usr/bin/Rscript --vanilla

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

Samples <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)
filepaths <- snakemake@input[["bam"]]

names(filepaths) <- gsub(".filtered.sortedByCoord.out.bam","",basename(filepaths))
filepaths <- filepaths[order(match(names(filepaths),row.names(Samples)))]

bamfiles <- BamFileList(filepaths)

gtffile <- snakemake@input[["gtf"]]
txdb <- makeTxDbFromGFF(gtffile,format="gtf",circ_seqs=character())
#save(txdb,file="txdb.RData")

end = as.logical(snakemake@params[["end"]])
frag = !end

ebg <- exonsBy(txdb,by="gene")

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=end,
                        ignore.strand=TRUE,
                        fragments=frag )

colData(se) <- DataFrame(Samples)
# dim(se)
# assayNames(se)
# head(assay(se),10)
# colSums(assay(se))
#For testing only for now
#save(se,file=args[4])

dds = DESeqDataSet(se,design=as.formula(snakemake@params[["model"]]))

dds <- dds[ rowSums(counts(dds)) > snakemake@params[["filt"]], ]

dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])
norm_counts = counts(dds, normalized=T)
write.table(data.frame("gene"=rownames(norm_counts), norm_counts), file=snakemake@output[[2]], sep='\t', row.names=F)