#This script installs the necessary packages to run ATACseqQC that conda can't install.
#To avoid installing in every ATACseqQC run (may cause conflicts)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
library(BiocManager)
install(snakemake@params[["BSgenome"]], update = TRUE, ask = FALSE)
install(snakemake@params[["Txdb"]], update = TRUE, ask = FALSE)
if (snakemake@params[["BSgenome"]] == "BSgenome.Mmusculus.UCSC.mm10") {
     library(GenomicScores)
     #This is to avoid an error over number of downloads in Rscript
     library(AnnotationHub)
     setAnnotationHubOption("MAX_DOWNLOADS", 100)
     gsco <- getGScores("phastCons60way.UCSC.mm10")
} else if (snakemake@params[["BSgenome"]] == "BSgenome.Hsapiens.UCSC.hg38") {
     install("phastCons100way.UCSC.hg38", update = TRUE, ask = FALSE)
} else if (snakemake@params[["BSgenome"]] == "BSgenome.Hsapiens.UCSC.hg19") {
     install("phastCons100way.UCSC.hg19", update = TRUE, ask = FALSE)
}