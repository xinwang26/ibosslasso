#######################################################################
#Author: Xin Wang
#Initial Date: Thu Jun  9 13:16:10 2016
#Program Description: this program generate and store the data for regression to save some time
#Completed: 
#Dependencies: 
#Main Reference: 
#Input: 
#Output:
#######################################################################
#######################################################################
#Function Description: data_gen generate a list of a data.table and one vector -- x,beta
#Dependencies: data.table
#Input: 3 random seed, scale of true beta
#Output: x, beta
#######################################################################
data_gen_cor_r <- function(seed1, n_full=n, p_full=p, rho=0.3){
  set.seed(seed1)
  sigma<-rbind(c(1,rho), c(rho,1))
  system.time(x <- as.data.table(matrix(c(mvrnorm(n=n_full*p_full/2, mu=c(0,0), Sigma=sigma)),ncol=p_full)))
  return(x)
}

data_gen_n<-function(seed1,n_full = n,p_full = p){
  #x
  set.seed(seed1)
  system.time(x <- as.data.table(matrix(rnorm(n_full*p_full),ncol=p_full)))
  x[,id:=.I]
  return(x)
}

data_gen_t <-function(seed1,n_full = n,p_full = p,df=2){
  #x
  set.seed(seed1)
  system.time(x <- as.data.table(matrix(rt(n_full*p_full,df),ncol=p_full)))
  x[,id:=.I]
  return(x)
}

data_gen_l <-function(seed1,n_full = n,p_full = p){
  #x
  set.seed(seed1)
  system.time(x <- as.data.table(matrix(rlnorm(n_full*p_full),ncol=p_full)))
  x[,id:=.I]
  return(x)
}

data_gen_m <-function(seed1,n_full = n,p_full = p){
  #x
  set.seed(seed1)
  system.time(x <- as.data.table(matrix(0.25*rlnorm(n_full*p_full)
                                        +0.25*c
                                        +0.25*rnorm(n_full*p_full)
                                        +0.25*rt(n_full*p_full,3)
                                        ,ncol=p_full)
                                 )
              )
  x[,id:=.I]
  return(x)
}

data_gen_vari<-function(seed1,n_full = n,p_full = p){
  #x
  set.seed(seed1)
  system.time(x <- as.data.table(matrix(rnorm(n_full*p_full)+runif(n_full*p_full,-10,10),ncol=p_full)))
  x[,id:=.I]
  return(x)
}

para_gen<-function(seed2,scale_b = scale_beta,p_full = p){
  #beta
  set.seed(seed2)
  beta_ind = sample(1:p_full,p_true)
  (beta_val = scale_beta + rnorm(p_true,0,scale_beta/5))
  beta = rep(0,p_full)
  beta[beta_ind] = beta_val
  beta = matrix(beta,length(beta))
  return(beta)
}


#cat(sum( sapply(ls(),function(x){object.size(get(x))})),"Bytes after random generate\n")
#y = ff(3*z[,1]+5*z[,15]+rnorm(n,0,1))


#system.time(x1 <- as.matrix(x))
#set.seed(seed3)
#  y = x1%*%beta+rnorm(n)
#write.csv(y,"y.csv",row.names = F)

data_gen_logistic<-function(seed1,n_full = n,p_full = p){
  
}