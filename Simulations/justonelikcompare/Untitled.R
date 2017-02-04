rm(list=ls())
library('mash')
set.seed(123)
c=chat_sim_fsingle_fixedomega(n = 10000,d = 8,omega = 2,esd = 0.1)
saveRDS(c,"chatfixedomega.rds")
diag(c$covmat)
source("~/matrix_ash/R/mashnobasescripts.R")
source("~/matrix_ash/R/truthscripts.R")
rm(A)

library("ExtremeDeconvolution")

R=8

c=readRDS("chatfixedomega.rds")
t=c$t;b=c$chat;se=c$shat
#hist(t)
#hist(b)
R=ncol(t)
L=diag(R)-1/R*as.vector(rep(1,R))%*%t(as.vector(rep(1,R)))
s.j=se/se

absmat=abs(t(L%*%t(c$t)))
index=which(rowSums(absmat>4)>0)
length(index)
mean(index<1000)

#index=which(rowMeans(abs(t(L%*%t(c$t))))>1.3)##choose the associations with average deviation large.

##show deviations
plot(rowSums(absmat>4))
points(index,rowSums(absmat>4)[index],col="blue")

sparsestrongt=t(L%*%t(t[index,]))

strongprojectedt=sparsestrongt


lvllist=genlvllist(s.j[index,],L = L[-1,])


##w should represent transformed betahats using all R subgroups and j genes
wfull=t(L%*%t(b))
sjmat=convertstandarderrors(s.j = se,L=L)##make this be the standard errors of all, it is the sqrt of the diagonal of LVL' for each j. Recall here L is still RxR so we can have the full set of deviations and thier standard errors to choose grid
dim(wfull)
dim(sjmat)

t.stat =sparsestrongt;v.j = lvllist;L = L[-1,];w = strongprojectedt[,-1]
init.cov.list=list()
init.cov.list[[1]]=cov(strongprojectedt)
#head(init.cov.list)
mean.mat=matrix(rep(0,R*R),ncol=R,nrow=R)  
ydata=w
xamp= 1
xcovar= init.cov.list
fixmean= TRUE     
ycovar=  v.j     
xmean=   mean.mat   
projection=list();for(l in 1:nrow(t.stat)){projection[[l]]=L}
e=extreme_deconvolution(ydata=ydata,ycovar=ycovar,xamp=xamp,xmean=xmean,xcovar=init.cov.list,fixmean=T,projection=projection)
se.train=se

covlist=e$xcovar
R=8;L=diag(R)-1/R*as.vector(rep(1,R))%*%t(as.vector(rep(1,R)))
train.b=data.frame(wfull)[,-1];se.train = se;covmat = covlist;A = A;pen = 1;L = L[-1,];J=nrow(train.b)
lik.mat=t(sapply(seq(1:J),function(j){
  V.gp.hat=diag(se.train[j,])^2;
  LVL=L%*%V.gp.hat%*%t(L);
  LSigL=L%*%covlist[[1]]%*%t(L)
  b.mle=train.b[j,]
  dmvnorm(x=b.mle, sigma=(LSigL + LVL),log=TRUE)}))

sum(log(exp(lik.mat)))
j=1
lik.mat[1]
V.gp.hat=diag(se.train[j,])^2;
LVL=L%*%V.gp.hat%*%t(L)
LSigL=L%*%covlist[[1]]%*%t(L)
dmvnorm(x=train.b[j,], sigma=(LSigL + LVL),log=TRUE)

```


