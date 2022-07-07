remove(list=ls())
library(lars)
library(glmnet)

library(ff)
library(ffbase)
library(biglars)

library(bigmemory)
library(bigalgebra)
library(biglasso)
library(irlba)

library(beepr)
library(mvtnorm)
n = 500
p = 1e3
n_true = 25
round =  1
#system.time(z <- (matrix(rnorm(1e6*1e2),ncol = 1e2)))

#system.time(md4<-biglars.fit(z,y,type = "lasso"))

#system.time(md5<-biglars.fit(z,y,type = "lasso"))
cat(sum(sapply(ls(),object.size)),"BYTES")
sstime=c()
true = as.ffdf(as.ff(rep(-100,p),dim=c(1,p)))
i=1
for(i in 1:round)
{1
  cat(sum(sapply(ls(),object.size)),"BYTES")
  #z <-as.ff(x = matrix(rnorm(n*p),n,p) ) 
  system.time(z <-(matrix(rnorm(n*p),ncol=p)))
  print("finished random X")
  cat(sum( sapply(ls(),function(x){object.size(get(x))})),"Bytes after random generate\n")
  #y = ff(3*z[,1]+5*z[,15]+rnorm(n,0,1))
  beta_ind = sample(1:p,n_true)
  beta_val = rgamma(n_true,5,1)
  beta = rep(0,p)
  beta[beta_ind] = beta_val
  beta = as.ff(
              matrix(rnorm(p),p,1) )
  #true = ffdfappend(true,beta,adjustvmode = F)
  # y = z%*%beta
  y = ffmatrixmult(z,beta)
  #sstime=c(system.time(assign(paste0("md",i),cv.glmnet(z,y,alpha=1)))[3],sstime)
  sstime=c(system.time(assign(paste0("md",i),biglars.fit(z,y,type = "lar")))[3],sstime)
  assign(paste0("md",i),coef(get(paste0("md",i)))[dim(coef(get(paste0("md",i))))[1],])
  cat(sum( sapply(ls(),function(x){object.size(get(x))})),"Bytes after fit\n")
  remove(z,y)
  gc()
  cat(sum( sapply(ls(),function(x){object.size(get(x))})),"Bytes after clear\n")
}
beep()
print(sstime)



#small data = iboss(big data)

#lasso(small data)


#result...?


as.bit.ff(ff(c(T,F,F)))
A =as.ff.bit((rep(FALSE,1000)))

