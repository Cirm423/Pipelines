log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
library(BiocManager)
install(snakemake@params[["BSgenome"]])
install(snakemake@params[["Txdb"]])
library(ATACseqQC)
bamfile <- snakemake@input[[1]]
bamfile.labels <- gsub(".bam", "", basename(bamfile))
#Change the outputs from number to names as you make the rule
pdf(snakemake@output[["fragmentSizeDistribution"]]) #fragmentSizeDistribution.pdf
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()
library(snakemake@params[["Txdb"]])
txs <- transcripts(snakemake@params[["Txdb"]])
library(snakemake@params[["BSgenome"]])
seqlev <- paste0("chr", c(1:21, "X", "Y"))
which <- as(seqinfo(snakemake@params[["BSgenome"]])[seqlev], "GRanges")
gal <- readBamFile(bamfile, which=which, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal) #outbam="shifted.bam" was there
pt <- PTscore(gal1, txs)
pdf(snakemake@output[["PTscore"]]) #"PTscore.pdf"
plot(mcols(pt)[, "log2meanCoverage"], mcols(pt)[, "PT_score"], 
    xlab="log2 mean coverage",
    ylab="Promoter vs Transcript")
dev.off()
nfr <- NFRscore(gal1, txs)
pdf(snakemake@output[["NFRscore"]]) #NFRscore.pdf
plot(mcols(nfr)[, "log2meanCoverage"], mcols(nfr)[, "NFR_score"], 
    xlab="log2 mean coverage",
    ylab="Nucleosome Free Regions score",
    main="NFRscore for 200bp flanking TSSs",
    xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()
tsse <- TSSEscore(gal1, txs)
pdf(snakemake@output[["TSSEscore"]]) #TSSEscore.pdf
hist(mcols(tsse)[, "TSS.enrichment.score"], breaks=100, 
main="Transcription Start Site (TSS) Enrichment Score", 
xlab="TSS enrichment score")
dev.off()
gc(reset=TRUE)
genome <- snakemake@params[["BSgenome"]]
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = ".")
rm(gal1)
gc(reset=TRUE)
library(ChIPpeakAnno)
bamfiles <- file.path(".",
                    c("NucleosomeFree.bam",
                    "mononucleosome.bam",
                    "dinucleosome.bam",
                    "trinucleosome.bam"))
pdf(snakemake@output[["cumulativePercentage"]]) #cumulativePercentage.pdf
cumulativePercentage(bamfiles[1:2], as(seqinfo(snakemake@params[["BSgenome"]])[seqlev], "GRanges"))
dev.off()
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
librarySize <- estLibSize(bamfiles)
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                    "mononucleosome",
                                    "dinucleosome",
                                    "trinucleosome")], 
                        TSS=TSS,
                        librarySize=librarySize,
                        seqlev=seqlev,
                        TSS.filter=0.5,
                        n.tile = NTILE,
                        upstream = ups,
                        downstream = dws)
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
    pdf(snakemake@output[["featureAligndHeatmap"]]) #featureAligndHeatmap.pdf
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                    zeroAt=.5, n.tile=NTILE)
    dev.off()
    out <- featureAlignedDistribution(sigs, 
                                reCenterPeaks(TSS, width=ups+dws),
                                zeroAt=.5, n.tile=NTILE, type="l", 
                                ylab="Averaged coverage")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
pdf(snakemake@output[["TSS_profile"]]) #TSS.profile.pdf
matplot(out, type="l", xaxt="n", 
    xlab="Position (bp)", 
    ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
    labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
library(MotifDb)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
pdf(snakemake@output[["CTCF_footprint"]]) #CTCF.footprint.pdf
sigs <- factorFootprints("shifted.bam", pfm=CTCF[[1]], 
                        genome=genome,
                        min.score="90%", seqlev=seqlev,
                        upstream=100, downstream=100)
dev.off()
pdf(snakemake@output[["CTCF_Vplot"]]) #CTCF.Vplot.pdf
vp <- vPlot("shifted.bam", pfm=CTCF[[1]], 
        genome=genome, min.score="90%", seqlev=seqlev,
        upstream=200, downstream=200, 
        ylim=c(30, 250), bandwidth=c(2, 1))
dev.off()
unlink("Rplots.pdf")