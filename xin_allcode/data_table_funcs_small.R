#######################################################################
#Author: Xin Wang
#Initial Date: Thu Jun  9 16:03:31 2016
#Program Description: This program define sorting functions in form of data.table
#Completed: 
#Dependencies: 
#Main Reference: 
#Input: 
#Output:
#######################################################################
library(data.table)
library(MASS)
#######################################################################
#Function Description: this function returns id of informative subset by dimension
#Dependencies: data.table
#Input: a data table
#Output: a vector of index
#######################################################################
all_vars_max <- function(x,y,sample_size=subsize,...)
{
  x1=x
  ni = ifelse("id" %in% names(x1),
              floor(sample_size/(dim(x1)[2]-1)), 
              floor(sample_size/(dim(x1)[2])))
  x1[,id:=.I]
  if(ni<1){stop("Wrong sample size!")}
  
  #need to deal with V1 separately
  ptm = proc.time()
  var1 = names(x1)[1]
  m_p =ifelse(sum(x1[,var1]>0)<dim(x1)[1]/2, ceiling(ni/2), floor(ni/2))
  #system.time(idx_p <- x1[,.I[order(-var1)][1:m_p],with=F])
  #system.time(idx_n <- x1[,.I[order(var1)][1:(ni-m_p)]])
  system.time(idx_p <- eval(parse(text=paste0("x1[,id[order(-",var1,")][1:",m_p,"]]"))))
  system.time(idx_n <- eval(parse(text=paste0("x1[,id[order(",var1,")][1:",ni-m_p,"]]"))))
  id_samp = c(idx_p,idx_n)
  proc.time()-ptm
  
  for(var in names(x1)[!names(x1)%in%c("id",var1)] ) #except V1 and id
  {
    ptm = proc.time()
    m_p =ifelse(sum(x1[,var,with=F]>0)<dim(x1)[1]/2, ceiling(ni/2), floor(ni/2))
    system.time(idx_p <- eval(parse(text=paste0("x1[-c(",paste(id_samp,collapse = ","),"),id[order(-",var,")][1:",m_p,"]]"))))
    system.time(idx_n <- eval(parse(text=paste0("x1[-c(",paste(id_samp,collapse = ","),"),id[order(",var,")][1:",ni-m_p,"]]"))))
    idx = c(idx_p,idx_n)
    id_samp = c(id_samp,idx)
    proc.time()-ptm
  }
  proc.time()-ptm
  return(id_samp)
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
  system.time(cor_resp <- unlist(lapply(x,cor,y=y)))
  system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
  sub_table =x[,to_select,with = F] 
  names(sub_table) = paste0("V",1:to_select_size)
  (selected = all_vars_max(sub_table))
}

#######################################################################
#Function Description: this function returns the observations have small value in un-correlated variables
#Dependencies: 
#Input: 
#Output:
#######################################################################
cor_var_min <-function(x,y,sample_size=subsize,prop = 0.5,...)
{
  (to_select_size = floor(dim(x)[2]*prop))
  system.time(cor_resp <- unlist(lapply(x[,names(x)[!names(x)%in%"id"],with=F],cor,y=y)))
  system.time(to_select <- sort(abs(cor_resp),index.return = T)$ix[1:to_select_size])
  sub_table =x[,to_select,with = F] 
  names(sub_table) = paste0("V",1:to_select_size)
  system.time(selected <- all_vars_max(sub_table))
  return(selected)
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
lev_sample <- function(x,y,sample_size = subsize,seed6 = seed6,...)
{
  #to compute leverage:
  X = cbind(matrix(rep(1,dim(x)[1]),ncol=1),x_mat)
  xx_inv <<- ginv(t(X)%*%X)
  cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F])
  x_lev = apply(cbind(v0=rep(1,dim(x)[1]),x[,-which(names(x)=="id"),with=F]),1,exact_lev)
  set.seed(seed6)
  idx = sample(1:dim(x)[1],size=sample_size,prob = x_lev)
  return(idx)
}



