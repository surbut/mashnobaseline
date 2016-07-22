
##process john's data'###

library("limma")
library("edgeR")

full <- read.table("~/Dropbox/john/counts_per_sample.txt", header = TRUE,
                   stringsAsFactors = FALSE)
full <- full[order(full$dir), ]
rownames(full) <- paste(full$dir, full$ind, full$bact, full$time, sep = ".")
counts <- t(full[, grep("ENSG", colnames(full))])
# Filter lowly expressed genes
counts <- counts[rowSums(cpm(counts) > 1) > 6, ]
dim(counts)
##now counts is JxN, transform to NxJ

counts=t(counts)

rep.col<-function(x,n){

  matrix(rep(x,each=n), ncol=n, byrow=TRUE)

}

voom2 <- function(counts){

  libsize.mat <- rep.col(rowSums(counts), dim(counts)[2]);

  voom.out <- log((counts+0.5), base=2) - log((libsize.mat+1), base=2)+ 6* log(10, base=2);

  return(voom.out)

}



groups <- full[, c("ind", "bact", "time", "extr", "rin")]
groups$bact <- gsub("\\+", "plus", groups$bact)
groups$ind <- factor(groups$ind)
groups$bact <- factor(groups$bact, levels = c("none", "Rv", "Rvplus", "GC",
                                              "BCG", "Smeg", "Yers", "Salm",
                                              "Staph"))
groups$time <- factor(groups$time, levels = c(4, 18, 48))
groups$extr <- factor(groups$extr)
head(groups)
dim(groups)
knownsamples=which(groups$time==18&groups$bact!="none")
counts_class=counts[knownsamples,]
class_labs=groups[knownsamples,"bact"]
voom_class <- voom2(counts_class);
model_mat <- model.matrix(~as.factor(class_labs)-1)[,-1] 
colSums(model_mat)
beta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2]);

sebeta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2])


for(k in 1:dim(model_mat)[2]){

      model_mat_temp <- cbind(model_mat[,k]);

      limma.obj <- limma::lmFit(t(voom_class), model_mat_temp,

                                weights=t(limma::voom(counts_class)$weights))





      beta_class[,k] <- as.matrix(limma.obj$coefficients[,1]);

      sebeta_class[,k] <- limma.obj$sigma*(as.matrix(limma.obj$stdev.unscaled[,1]));
      

    }


