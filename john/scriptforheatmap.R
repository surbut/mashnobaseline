batchcorrectdata=read.table("~/Dropbox/oldnobaseline/john/s1table.txt",header = T,stringsAsFactors = F)[,-c(1,2)]
anno=read.table("~/Dropbox/oldnobaseline/john/annotation.txt",header=T)
library('limma')

library('gplots')
library('colorRamps')

groups <- anno[, c("ind", "bact", "time", "extr", "rin")]
groups$bact <- gsub("\\+", "plus", groups$bact)
groups$ind <- factor(groups$ind)
groups$bact <- factor(groups$bact, levels = c("none", "Rv", "Rvplus", "GC",
"BCG", "Smeg", "Yers", "Salm",
"Staph"))
groups$time <- factor(groups$time, levels = c(4, 18, 48))
groups$extr <- factor(groups$extr)
head(groups)

mash.means=read.table("batchcorrectwithEDposterior.means.txt")[,-1]
colnames(mash.means)=levels(groups$bact)[-1]

library('gplots')
library('colorRamps')
x=as.matrix(mash.means)

png("estimatedeffectsJohn.png")
heatmap.2(x/max(x),
dendrogram = "none",col=blue2green,density.info = "none",trace="none",main=paste0("EstimatedEffects, MASH"))
dev.off()


x=as.matrix(mash.means)

png("estimatedeffectsJohn.png")
heatmap.2(x/max(diag(x)),
#Rowv = NULL,Colv = NULL,
dendrogram = "none",col=blue2green,density.info = "none",trace="none",main=paste0("EstimatedEffects, MASH"))
dev.off()
