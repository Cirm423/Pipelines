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
            covariates = covariates,
            dbtype = "tabix",
            dbdir = snakemake@output[["db"]]
            )
} else {
    methDB=methRead(files.list,
            sample.id=as.list(row.names(samples)),
            assembly=snakemake@params[["assembly"]],
            treatment=snakemake@params[["treatment"]],
            context="CpG",
            mincov = snakemake@params[["mincov"]],
            pipeline = pipeline,
            dbtype = "tabix",
            dbdir = snakemake@output[["db"]]
            )
}

#The following 2 plots are for single samples, so iterate over them
dir.create(snakemake@output[["CpG_methylation"]], recursive = TRUE)
dir.create(snakemake@output[["CpG_coverage"]], recursive = TRUE)

for (sample_index in 1:dim(samples)[1]) {
    sample_name <- row.names(samples)[sample_index]
    #Print and save CpG methylation % plot
    out_file_meth <- paste(snakemake@output[["CpG_methylation"]],"/",sample_name,"-methylation.pdf",sep="")
    pdf(file=out_file_meth)
    getMethylationStats(methDB[[sample_index]],plot=TRUE,both.strands=FALSE)
    dev.off()
    #Print and save CpG coverage plot
    out_file_cov <- paste(snakemake@output[["CpG_coverage"]],"/",sample_name,"-coverage.pdf",sep="")
    pdf(file=out_file_cov)
    getCoverageStats(methDB[[sample_index]],plot=TRUE,both.strands=FALSE)
    dev.off()
}


#Save RDS for manual use
saveRDS(methDB, file=snakemake@output[["RDS"]])

#Filter samples based on read coverage (might be useful to change this to an option)
filtered.methDB=filterByCoverage(methDB,lo.count=snakemake@params[["mincov"]],lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)


#Merge samples for further downstream analysis (destrand is only useful in CpG, which is the only case for now)
if (snakemake@params[["window_s"]] > 1) {
    tiles = tileMethylCounts(filtered.methDB,win.size=snakemake@params[["window_s"]],step.size=snakemake@params[["step_s"]],cov.bases = snakemake@params[["tile_cov"]],mc.cores=snakemake@threads[[1]])
    meth = unite(tiles, destrand=FALSE, min.per.group=snakemake@params[["min_group"]], mc.cores=snakemake@threads[[1]])
} else {
    meth=unite(filtered.methDB, destrand=FALSE, min.per.group=snakemake@params[["min_group"]], mc.cores=snakemake@threads[[1]])
}

#Print and save correlation plot
pdf(file=snakemake@output[["correlation"]])
getCorrelation(meth,plot=TRUE)
dev.off()

#Print and save clustering plot
pdf(file=snakemake@output[["cluster"]])
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#Print and save PCA screen plot
pdf(file=snakemake@output[["screen"]])
PCASamples(meth, screeplot=TRUE)
dev.off()

#Print and save PCA plot
pdf(file=snakemake@output[["PCA"]])
PCASamples(meth)
dev.off()

#Differential analysis with overdispersion, depending on the presence of covariates
if (length(colnames(samples)) > 1) {
    myDiff<-calculateDiffMeth(meth,
                                    covariates=covariates,
                                    overdispersion="MN",mc.cores=snakemake@threads[[1]])
} else {
    myDiff<-calculateDiffMeth(meth,
                                    overdispersion="MN",mc.cores=snakemake@threads[[1]])
}

# Report methylation and save them
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
write.table(getData(myDiff25p.hyper), file=snakemake@output[["hyper"]], sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
write.table(getData(myDiff25p.hypo), file=snakemake@output[["hypo"]], sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
write.table(getData(myDiff25p), file=snakemake@output[["all_diff"]], sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

chrDiff25p=diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
write.table(getData(chrDiff25p), file=snakemake@output[["chr_diff"]], sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# Create bedgraph of all diff, don't save directly to file to avoid the track line for bigWig conversion
#Disable scientific notation in bedgraph output.
options(scipen=999)

bg = bedgraph(myDiff, col.name = "meth.diff", unmeth = FALSE)
write.table(bg, file=snakemake@output[["bed"]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#Annotation of differentially methylated stuff
gene.obj=readTranscriptFeatures(snakemake@input[["annot"]],remove.unusual=FALSE)

anno = annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

#Create the annotation dataframe
anno.df = merge(getAssociationWithTSS(anno),getMembers(anno),by.x=c("target.row"), by.y=0)

write.table(anno.df, file=snakemake@output[["annotation"]], sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)