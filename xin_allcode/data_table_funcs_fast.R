#######################################################################
#Author: Xin Wang
#Initial Date: Thu Jun  9 16:03:31 2016
#Program Description: This program define sorting functions in form of data.table
#Completed: 
#Dependencies: sortcpp.R(contains cpp function getidx)
#Main Reference: 
#Input: 
#Output:
#######################################################################
library(data.table)
library(MASS)
library(Rcpp)
library(inline)
source("sortcpp.R",echo =T,max.deparse.length = 10000)
#######################################################################
#Function Description: this function returns id of informative subset by dimension
#Dependencies: data.table
#Input: a data table
#Output: a vector of index
#######################################################################
all_vars_max_real <- function(x,y=y,sample_size=subsize,...)
{
  if(dim(x)[1]<=sample_size){return(x[,id])}
  ni = ifelse("id" %in% names(x),
              floor(sample_size/(dim(x)[2]-1)/2), #the fastsort function will take the half
              floor(sample_size/(dim(x)[2])/2))
  x = data.table(x)
  x[,id:=.I]
  if(ni<1){stop("Wrong sample size!")}
  ######## unit of id get
  idx <- c()
  idx <- x[order(x[,V1])[c(1:ni,(dim(x)[1]-ni+1):dim(x)[1])],id]
  count_unique =0
  #if(length(idx)>2*ni){idx = idx[c(1:ni,(length(idx)-ni+1):length(idx))]}
  for(j in 2:(dim(x)[2]-1)) {#subtract the 1 dim for id
    ## tmp <- .Call("get", Z[-idx.oD,j], r)
    aaa =eval(parse(text=paste0("x[-idx,.(V",j,",id)]")))
    tmp <-  order(aaa[,1])[c(1:ni,(dim(aaa)[1]-ni+1):dim(aaa)[1])]
    tmp = aaa[tmp,id]
    #cat(tmp)
    #if(length(tmp)>1*ni){tmp = tmp[c(1:(2*ni))]}
    idx1 <- c(idx, tmp)
    count_unique = count_unique + (length(idx1)==length(unique(idx1)))
    idx = idx1
    #cat(idx)
  }
  return(idx)
}
cor_var_max_real <-function(x,y,sample_size=subsize,prop = 0.2,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  if (sum(apply(x,2,var)==0)>0){
    corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-(which(apply(x,2,var)==0)),with=F],cor,y=y)) )
  }
  else{
    corcalc_time <- system.time(cor_resp <- unlist(lapply(x,cor,y=y)) )
  }
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  sub_table = data.table(x[,c(to_select),with=F]) #then require x must contains id
  names(sub_table) = c(paste0("V",1:to_select_size))
  selected <- all_vars_max_real(x = sub_table,y=y, sample_size)
  return(selected)
}
#######################################################################
#Function Description: this function returns id of informative subset by dimension
#Dependencies: data.table
#Input: a data table
#Output: a vector of index
#######################################################################
all_vars_max <- function(x,y=y,sample_size=subsize,...)
{
  if(dim(x)[1]<=sample_size){return(x[,id])}
  ni = ifelse("id" %in% names(x),
              floor(sample_size/(dim(x)[2]-1)/2), #the fastsort function will take the half
              floor(sample_size/(dim(x)[2])/2))
  #x[,id:=.I]
  if(ni<1){stop("Wrong sample size!")}
  ######## unit of id get
  idx <- c()
  idx <- getIdx(x[,V1], 
                x[,V1], 
                ni, dim(x)[1])#x f r max
  for(j in 2:(dim(x)[2]-1)) {#subtract the 1 dim for id
    ## tmp <- .Call("get", Z[-idx.oD,j], r)
    tmp <- getIdx(eval(parse(text=paste0("x[-idx,V",j,"]"))), #everytime exclude selected and get new
                  eval(parse(text=paste0("x[,V",j,"]"))), 
                  ni,dim(x)[1])
    idx <- unique(c(idx, tmp))
  }
  return(idx)
}
#######################################################################
#Function Description: this function returns the observations having largest values in those variables most correlated to response
#Dependencies: 
#Input: 
#Output:
#######################################################################
cor_var_max <-function(x,y,sample_size=subsize,prop = 0.3,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  sub_table =x[,c(to_select,dim(x)[2]),with = F] #then require x must contains id
  names(sub_table) = c(paste0("V",1:to_select_size),"id")
  selected <- all_vars_max(sub_table)
  return(selected)
}
cor_var_max_sis <-function(x,y,sample_size=subsize,prop = 0.3,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  sub_table =x[,c(to_select,dim(x)[2]),with = F] #then require x must contains id
  names(sub_table) = c(paste0("V",1:to_select_size),"id")
  select_time = system.time(selected <- all_vars_max(sub_table))
  return(list(rows= selected,cols=to_select,select_time =select_time[3]))
}

#######################################################################
#Function Description: this function returns to the observations with largest total abslute value of correlation coefficient with the responses
#Dependencies: 
#Input: 
#Output:
#######################################################################
sum_prod<-function(corr_vec,x_vec){return(sum(corr_vec,x_vec))}
cor_to_max <-function(x,y,sample_size=subsize,...){
  all_pred = x[,names(x)[!names(x)%in%"id"],with=F]
  system.time(cor_resp <- unlist(lapply(all_pred,cor,y=y)))
  system.time(sum_cor <- apply(abs(all_pred),1,sum_prod,corr_vec=abs(cor_resp)))
  samp_id  <- sort(sum_cor,decreasing = T,index.return = T)$ix[1:sample_size]
  return(samp_id)
  }


#######################################################################
#Function Description: this function returns the observations having larges mods
#Dependencies: 
#Input: 
#Output:
#######################################################################
norms <- function(xvec){return(sum(xvec^2))}
mod_var_max <-function(x,y,sample_size,...)
{
  x1 = x[,names(x)[!names(x)%in%"id"],with=F]
  x1$norm = apply(x1,1,norms)
  ni = sample_size
  system.time(sample_id <- x1[,.I[order(-norm)][1:ni]])
  return(sample_id)
}

#######################################################################
#Function Description: this function calculate exact leveraging score for targeted row
#Dependencies: x_mat-- matrix/bigmatrix formatted design matrix
#Input: by global environment to save memory -- (X'X)^-1 will be xx_inv from loop
#Output: index of leveraging sample
#######################################################################

exact_lev <- function(x_row){
  x_row = as.numeric(x_row)
  return(x_row%*%xx_inv%*%x_row)
}

roundUp <- function(x) 10^ceiling(log10(x))
approx_lev <- function(X,p1=p)
{
  r1 <- 501
  r1 = min(r1,p1+1)
  (r2 <- round(floor(log(n)*4),-1) )
  SX <- matrix(0,r1,p1+1)
  rn <- rmultinom(1,n,prob=rep(1,r1))
  for(i in 1:r1){
    index.pm <- sample(1:n,rn[i],F)  #the randomly seleted actual data part index
    Dd <- 2*rbinom(rn[i],1,0.5)-1 #the diagonal D
    SX[i,] <- Dd%*%(matrix(X[index.pm,],,p1+1))  #the DH? --> the SRHT
  }
  sv.cw <- svd(SX); #SX equiv to Pi_1*A in the paper?
  V <- sv.cw$v
  D <- sv.cw$d
  Rinv <- V%*%diag(D);
  SRinv <- matrix(0,r2,p1+1);
  rn <- rmultinom(1,r1,prob=rep(1,r2));
  Rinv <- t(Rinv); #still use FLJT not LJT
  for(i in 1:r2){
    index.pm <- sample(1:r1,rn[i],F);
    Dd <- 2*rbinom(rn[i],1,0.5)-1;
    SRinv[i,] <- Dd%*%(matrix(Rinv[index.pm,],,r1));
  }
  B <- X%*%t(SRinv);
  p1.slev <- rep(0,n);
  for(i in 1:n) {
    p1.slev[i]  <-  sqrt(B[i,]%*%B[i,]);#the leverage score
  }
  p1.slev  <-  p1.slev/sum(p1.slev)
  return(p1.slev)
}

lev_sample <- function(x,y,sample_size = subsize,seed6 = seed6,...)
{
  #to compute leverage:
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat)
  xx_inv <<- ginv(t(X)%*%X)
  #cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F])
  x_lev = apply(cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F]),1,exact_lev)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  return(idx)
}
lev_sample_real <- function(x,y,sample_size = subsize,seed6 = seed6,...)
{
  #to compute leverage:
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat)
  xx_inv <<- ginv(t(X)%*%X)
  #cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F])
  x_lev = apply(cbind(v0=rep(1,dim(x)[1]),x_mat),1,exact_lev)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  return(idx)
}
approxlev_sample <- function(x,y,sample_size = subsize,seed6 = seed6,...)
{
  #to compute leverage:
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat)
  #cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F])
  x_lev = approx_lev(X)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  return(idx)
}

cor_lev_sample <-function(x,y,sample_size=subsize,prop = 0.3,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat[,to_select])
  xx_inv <<- ginv(t(X)%*%X)
  #cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F])
  x_lev = apply(X,1,exact_lev)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  return(idx)
}

cor_lev_sample_sis <-function(x,y,sample_size=subsize,prop = 0.3,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat[,to_select])
  tic = proc.time()
  xx_inv <<- ginv(t(X)%*%X)
  #cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F])
  x_lev = apply(X,1,exact_lev)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  select_time = proc.time() - tic
  return(list(rows=idx,cols=to_select,select_time=select_time[3]))
}



cor_alev_sample <-function(x,y,sample_size=subsize,prop = 0.3,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat[,to_select])
  x_lev =approx_lev(X,p1 = dim(X)[2]-1)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  return(idx)
}

cor_alev_sample_sis <-function(x,y,sample_size=subsize,prop = 0.3,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
  corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat[,to_select])
  tic = proc.time()
  x_lev =approx_lev(X,p1 = dim(X)[2]-1)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  select_time = proc.time() - tic
  return(list(rows=idx,cols=to_select,select_time=select_time[3]))
}
#######################################################################
#Function Description:  This function implement the data selection for logistic two step innovative approach
#Dependencies: 
#Input: whole data set through global env
#Output: the selected subset
#######################################################################

min_fun<-function(x,num_covariate=p_true)
{value=1/(x^2*(exp(x)/((1+exp(x))^2))^(num_covariate+1));return(value)}

two_step_range<-function(x
                         ,y=y
                         ,rand_size = r0
                         ,sample_size=subsize
                         ,...
  
)
{
  #first step
  set.seed(seed7)
  rsamp<- sample(1:dim(x)[1],r0) 
  rand_x <- as.matrix(x_mat[rsamp,])
  rand_y <-  as.matrix(y[rsamp])
  md0_cv <- cv.glmnet(rand_x,as.factor(rand_y),alpha=islasso,family = "binomial")
  (p_selected=sum(coef(md0_cv)!=0))
  op_result=optimize(min_fun,num_covariate=2*p_true,interval=c(-10,10))#currently there is no better solution 
  #it is a conflict that optimization need to know number of covariate while number of covariates need to be decided in the last step of lasso 
  #assuming no information about true variables available, it is reasonable to use selected variables from the first simple random sample
  c_star=op_result$minimum
  ci <- cbind(1,x_mat)%*%coef(md0_cv)
  delta= quantile(ci,0.6)- quantile(ci,0.2)
  ci =as.data.table(as.matrix(ci))
  names(ci)="xb"
  B = ci[(xb>c_star-delta)&(xb<c_star+delta),.I]
  subtable <- x[B,]
  system.time(selected <- all_vars_max(subtable))
  return(selected)
}
  




