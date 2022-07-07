#######################################################################
#Author: 
#Initial Date: Wed Aug 30 23:15:33 2017
#Program Description: Original Approximate Leveraging method 
#Completed: 
#Dependencies: 
#Main Reference: 
#Input: 
#Output:
#######################################################################

roundUp <- function(x) 10^ceiling(log10(x))
approx_lev <- function(X)
{
r1 <- p+1
(r2 <- round(floor(log(n)*4),-1) )
SX <- matrix(0,r1,p+1)
rn <- rmultinom(1,n,prob=rep(1,r1))
for(i in 1:r1){
  index.pm <- sample(1:n,rn[i],F)  #the randomly seleted actual data part index
  Dd <- 2*rbinom(rn[i],1,0.5)-1 #the diagonal D
  SX[i,] <- Dd%*%(matrix(X[index.pm,],,p+1))  #the DH? --> the SRHT
}
sv.cw <- svd(SX); #SX equiv to Pi_1*A in the paper?
V <- sv.cw$v
D <- sv.cw$d
Rinv <- V%*%diag(D);
SRinv <- matrix(0,r2,p+1);
rn <- rmultinom(1,p+1,prob=rep(1,r2));
Rinv <- t(Rinv); #still use FLJT not LJT
for(i in 1:r2){
  index.pm <- sample(1:(p+1),rn[i],F);
  Dd <- 2*rbinom(rn[i],1,0.5)-1;
  SRinv[i,] <- Dd%*%(matrix(Rinv[index.pm,],,p+1));
}
B <- X%*%t(SRinv);
p.slev <- rep(0,n);
for(i in 1:n) {
  p.slev[i]  <-  sqrt(B[i,]%*%B[i,]);#the leverage score
}
p.slev  <-  p.slev/sum(p.slev)
return(p.slev)
}

idx.lv  <-  sample(1:n, k, T, p.slev)


