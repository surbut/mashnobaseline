
library('clusterGeneration')
set.seed(123)

library("MASS")
library("mvtnorm")
n=1000;d=5;betasd=1;esd=0.1;K=10
#chatsimtest=function(n=1000,d=3,betasd=1,esd=0.1,K=10){
  J=0.10*n
print(J)
  z = sample(K,J,replace=TRUE)
print(z)# randomly sample factor to be loaded on for each real snp
#rm(.Random.seed)
print(.Random.seed)
  covmat=lapply(seq(1:K),function(k){
    A=genPositiveDefMat("eigen",dim=d)$Sigma
    A/max(diag(A))##scale so max diag is 1
  })
print(covmat)
  mus=rnorm(n)  ###generate a list of n mus
print(mus[1:10])
  mumat=matrix(rep(mus,d),ncol=d) ##generate a matrix of mus for each gene
beta=t(sapply(seq(1:J),function(j){
    k=z[j]
    omega=abs(rnorm(1,mean=0,sd=betasd))##effect size variance can be big or small
    print(omega)
    mvrnorm(1,mu=rep(0,d),Sigma=omega*covmat[[k]])
    
  }))
print(beta[1:5,1:5])



  print(runif(1))
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  print(runif(1))
  c=beta+mumat
  print(runif(1))
  sj=abs(matrix(rnorm(n*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
  t=chat/sj
  return(list(beta=beta,chat=chat,covmat=covmat,components=z,t=t,mumat=mumat,shat=sj,error=e,ceff=c))
}

set.seed(123)
f=chatsimtest(n = 1000,d = 5,betasd = 1,esd = 0.1)
f$beta[1,]

set.seed(1)
