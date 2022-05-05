log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
library(BiocManager)
install(snakemake@params[["BSgenome"]])
install(snakemake@params[["Txdb"]])
library(ATACseqQC)
bamfile <- snakemake@input[[1]]
bamfile.labels <- gsub(".bam", "", basename(bamfile))
pdf(NULL) #This supposedly stops Rplots.pdf from being created
pdf(snakemake@output[["fragmentSizeDistribution"]]) #fragmentSizeDistribution.pdf
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()
library(snakemake@params[["Txdb"]],character.only=TRUE)
snake_txdb = get(snakemake@params[["Txdb"]])
txs <- transcripts(snake_txdb)
library(snakemake@params[["BSgenome"]],character.only=TRUE)
snake_BS = get(snakemake@params[["BSgenome"]])
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
library(Rsamtools)
bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
if (snakemake@params[["BSgenome"]] == "BSgenome.Mmusculus.UCSC.mm10"){
     seqlev <- paste0("chr", c(1:19, "X", "Y"))
     library(GenomicScores)
     #This is to avoid an error over number of downloads in Rscript
     library(AnnotationHub)
     setAnnotationHubOption("MAX_DOWNLOADS", 100)
     gsco <- getGScores("phastCons60way.UCSC.mm10")
} else if (snakemake@params[["BSgenome"]] == "BSgenome.Hsapiens.UCSC.hg38"){
     seqlev <- paste0("chr", c(1:21, "X", "Y"))
     install("phastCons100way.UCSC.hg38")
     library(phastCons100way.UCSC.hg38)
     gsco <- phastCons100way.UCSC.hg38
} else if (snakemake@params[["BSgenome"]] == "BSgenome.Hsapiens.UCSC.hg19"){
     seqlev <- paste0("chr", c(1:21, "X", "Y"))
     install("phastCons100way.UCSC.hg19")
     library(phastCons100way.UCSC.hg19)
     gsco <- phastCons100way.UCSC.hg19
}
which <- as(seqinfo(snake_BS)[seqlev], "GRanges")
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
dir.create(snakemake@params[["path"]])
shiftedBamfile <- file.path(snakemake@params[["path"]], "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
pt <- PTscore(gal1, txs)
pdf(snakemake@output[["PTscore"]]) #"PTscore.pdf"
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()
nfr <- NFRscore(gal1, txs)
pdf(snakemake@output[["NFRscore"]]) #NFRscore.pdf
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()
tsse <- TSSEscore(gal1, txs)
pdf(snakemake@output[["TSSEscore"]]) #TSSEscore.pdf
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()
gc(reset=TRUE)
genome <- snake_BS
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, conservation=gsco, outPath = snakemake@params[["path"]])
rm(gal1)
gc(reset=TRUE)
library(ChIPpeakAnno)
bamfiles <- file.path(snakemake@params[["path"]],
                    c("NucleosomeFree.bam",
                    "mononucleosome.bam",
                    "dinucleosome.bam",
                    "trinucleosome.bam"))
pdf(snakemake@output[["cumulativePercentage"]]) #cumulativePercentage.pdf
cumulativePercentage(bamfiles[1:2], as(seqinfo(snake_BS)[seqlev], "GRanges"))
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
sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                        genome=genome,
                        min.score="90%", seqlev=seqlev,
                        upstream=100, downstream=100)
dev.off()
pdf(snakemake@output[["CTCF_footprint"]]) #CTCF.footprint.pdf
sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                        genome=genome,
                        min.score="90%", seqlev=seqlev,
                        upstream=100, downstream=100)
dev.off()
pdf(snakemake@output[["CTCF_Vplot"]]) #CTCF.Vplot.pdf
vp <- vPlot(shiftedBamfile, pfm=CTCF[[1]], 
        genome=genome, min.score="90%", seqlev=seqlev,
        upstream=200, downstream=200, 
        ylim=c(30, 250), bandwidth=c(2, 1))
dev.off()
#unlink("Rplots.pdf")