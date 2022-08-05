log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(methylKit)
library(genomation)

pdf(NULL) #This supposedly stops Rplots.pdf from being created

samples <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

files <- snakemake@input[["meth"]]

#Adapting to the different files between bismark and methydackel
if (snakemake@params[["mode"]] == "bismark") {
    names(files) <- sapply(strsplit(basename(files), split="-pe|-se"), "[[", 1)
    #Making sure the files are in the same order as the samples + treatment vector
    files <- files[order(match(names(files), row.names(samples)))]
    pipeline <- "bismarkCoverage"
} else {
    #Change the split pattern for the correct methyldackel one.
    names(files) <- sapply(strsplit(basename(files), "_CpG.methylKit", fixed = TRUE), "[[", 1)
    #Making sure the files are in the same order as the samples + treatment vector
    files <- files[order(match(names(files), row.names(samples)))]
    pipeline <- "amp"
}

#Convert files to list as needed for the functions

files.list <- as.list(files)

#Read the files into a db, depending on covariates or not
if (length(colnames(samples)) > 1) {
    covariates = samples[,-1]
    methDB=methRead(files.list,
            sample.id=as.list(row.names(samples)),
            assembly=snakemake@params[["assembly"]],
            treatment=snakemake@params[["treatment"]],
            context="CpG",
            mincov = snakemake@params[["mincov"]],
            pipeline = pipeline,
            covariates = covariates
            )
} else {
    methDB=methRead(files.list,
            sample.id=as.list(row.names(samples)),
            assembly=snakemake@params[["assembly"]],
            treatment=snakemake@params[["treatment"]],
            context="CpG",
            mincov = snakemake@params[["mincov"]],
            pipeline = pipeline
            )
}

#Filter samples based on read coverage (might be useful to change this to an option)
filtered.methDB=filterByCoverage(methDB,lo.count=snakemake@params[["mincov"]],lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

#Generate bedgraphs depending on tile analysis or sCpG
if (snakemake@params[["window_s"]] > 1) {
    tiles = tileMethylCounts(filtered.methDB,win.size=snakemake@params[["window_s"]],step.size=snakemake@params[["step_s"]],cov.bases = snakemake@params[["tile_cov"]],mc.cores=snakemake@threads[[1]])
    bg = bedgraph(tiles,col.name="perc.meth", unmeth=FALSE)
    for (sample_index in 1:dim(samples)[1]) {
        write.table(bg[[sample_index]], file=snakemake@output[["bed"]][sample_index], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
} else {
    bg = bedgraph(filtered.methDB,col.name="perc.meth", unmeth=FALSE)
    for (sample_index in 1:dim(samples)[1]) {
        write.table(bg[[sample_index]], file=snakemake@output[["bed"]][sample_index], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
}