#! /usr/bin/Rscript --vanilla

## /home/share/dcyleung/Program/Pipeline/For_R/pois-test_RSEM.R

args = commandArgs(TRUE)

alt = args[1]

a <- read.table("results/single/counts.temp",sep="\t",header=FALSE,quote="")

b <- as.matrix(data.frame(ceiling(a$V5),ceiling(a$V6)))

#write.table(b,sep="\t",file="test",quote=FALSE,row.names=FALSE,col.names=FALSE)
x <- vector(mode="numeric", length=0)
for(row in 1:nrow(b)) {x <- append(x,poisson.test(as.vector(b[row,]),c(1,1),alternative=alt)$p.value)}

c <- data.frame(x)
a2 <- data.frame(a,V12=c$x)
a3 <- data.frame(a2,V13=p.adjust(a2$V12,method="fdr",n=length(a2$V12)))
write.table(a3,sep="\t",file="results/single/counts_diff.temp",quote=FALSE,row.names=FALSE,col.name=FALSE)


