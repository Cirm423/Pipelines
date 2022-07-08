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
} else {
    #Change the split pattern for the correct methyldackel one.
    names(files) <- sapply(strsplit(basename(files), "_CpG.methylKit", fixed = TRUE), "[[", 1)
    #Making sure the files are in the same order as the samples + treatment vector
    files <- files[order(match(names(files), row.names(samples)))]
}

#Read the files into a db
methDB=methRead(files,
           sample.id=as.list(row.names(samples)),
           assembly=snakemake@params[["assembly"]],
           treatment=snakemake@params[["treatment"]],
           context="CpG",
           mincov = snakemake@params[["mincov"]],
           dbtype = "tabix",
           dbdir = "results/methylDB"
           )

#Print and save CpG methylation % plot
pdf(file=snakemake@output[["CpG_methylation"]])
getMethylationStats(methDB,plot=TRUE,both.strands=FALSE)
dev.off()

#Print and save CpG coverage plot
pdf(file=snakemake@output[["CpG_coverage"]])
getCoverageStats(methDB,plot=TRUE,both.strands=FALSE)
dev.off()

#Filter samples based on read coverage (might be useful to change this to an option)
filtered.methDB=filterByCoverage(methDB,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

#Merge samples for further downstream analysis (destrand is only useful in CpG, which is the only case for now)
meth=unite(methDB, destrand=TRUE, min.per.group=snakemake@params[["min_group"]], mc.cores=snakemake@threads[[1]])

#Print and save correlation plot
pdf(file=snakemake@output[["correlation"]])
getCorrelation(meth,plot=TRUE)
dev.off()

#Print and save clustering plot
pdf(file=snakemake@output[["cluster"]])
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()

#Print and save PCA screen plot
pdf(file=snakemake@output[["screen"]])
PCASamples(meth, screeplot=TRUE)
dev.off()

#Print and save PCA plot
pdf(file=snakemake@output[["PCA"]])
PCASamples(meth)
dev.off()

#Checking if there are covariates in the samples table, if then do covariate analysis and actual differential analysis with overdispersion
if (length(colnames(samples)) > 1) {
    covariates <- samples[,-1]
    sim.methylBase<-dataSim(replicates=length(row.names(samples)),sites=1000,
                            treatment=snakemake@params[["treatment"]],
                            covariates=covariates,
                            sample.ids=row.names(samples)
                            )
    myDiff<-calculateDiffMeth(sim.methylBase,
                                    covariates=covariates,
                                    overdispersion="MN",mc.cores=snakemake@threads[[1]])
} else {
    sim.methylBase<-dataSim(replicates=6,sites=1000,
                            treatment=nakemake@params[["treatment"]],
                            sample.ids=row.names(samples)
                            )

    myDiff<-calculateDiffMeth(sim.methylBase,
                                    overdispersion="MN",mc.cores=snakemake@threads[[1]])
}

###Find a way to save these things into files

# Report methylation
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

#Annotation of differentially methylated stuff
gene.obj=readTranscriptFeatures(snakemake@input[["annot"]])

annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)