library(ExtremeDeconvolution)
library(mash)

c=readRDS("chatfixedomega.rds")
t=c$t;b=c$chat;se=c$shat;R=ncol(t)
L=diag(R)-1/R*as.vector(rep(1,R))%*%t(as.vector(rep(1,R)))
wfull=t(L%*%t(b))
s.j=se/se
absmat=abs(t(L%*%t(c$t)))
index=which(rowSums(absmat>4)>0)
factor.mat=as.matrix(read.table("sparseF_F.out"))
lambda.mat=as.matrix(read.table("sparseF_lambda.out"))
strongprojectedtsimulations=t(L%*%t(t[index,]))
t.stat = strongprojectedtsimulations;factor.mat = factor.mat;v.j = s.j;lambda.mat = lambda.mat;K = 3;P = 3

R=ncol(t.stat)

init.cov=init.covmat(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K=K,P=P)
init.cov.list=list()
for(i in 1:K){init.cov.list[[i]]=init.cov[i,,]}


init.cov2=init.covmat(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K=K,P=P)
init.cov.list2=list()
for(i in 1:K){init.cov.list2[[i]]=init.cov2[i,,]}
all.equal(init.cov.list,init.cov.list2)
mean.mat=matrix(rep(0,ncol(t.stat)*nrow(t.stat)),ncol=ncol(t.stat),nrow=nrow(t.stat))
ydata=  t.stat
xamp= rep(1/K,K)
xcovar= init.cov.list
fixmean= TRUE
ycovar=  v.j
xmean=   mean.mat
projection= list();for(l in 1:nrow(t.stat)){projection[[l]]=diag(1,R)}
e1=extreme_deconvolution(ydata=ydata,ycovar=ycovar,xamp=xamp,xmean=xmean,xcovar=init.cov.list,fixmean=T,projection=projection)
e2=extreme_deconvolution(ydata=ydata,ycovar=ycovar,xamp=xamp,xmean=xmean,xcovar=init.cov.list,fixmean=T,projection=projection)
cat("Identity check\n")
print(identical(e1,e2))
