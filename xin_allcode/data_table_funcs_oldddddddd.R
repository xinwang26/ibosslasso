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
#######################################################################
#Function Description: this function returns subset of rows of largest values for each variables in order of columns comes with data
#Dependencies: data.table
#Input: a data table
#Output: a vector of index
#######################################################################
all_vars_max <- function(x,sample_size=subsize)
{
  x1 = data.table(x)
  x1[,id:=.I]
  ptm = proc.time()
  ni = floor(sample_size/(dim(x1)[2]-1))
  if(ni<1){stop("Wrong sample size!")}
  samp_info = c()
  id = c()
  #NEED to deal with the first dimension -- V1
  m_p =ifelse(sum(x[,var,with=F]>0)<dim(x)[1]/2, ceiling(ni/2), floor(ni/2))
  system.time(idx_p <- x1[,.I[order(-V1)][1:m_p]])
  system.time(idx_n <- x1[,.I[order(V1)][1:(ni-m_p)]])
  id = c(idx_p,idx_n)
  samp_info = rbind(samp_info,x1[id,])
  for(var in names(x1)[!names(x1)%in%c("V1","id")] ) #except V1 and id
  {
    m_p =ifelse(sum(x[,var,with=F]>0)<dim(x)[1]/2, ceiling(ni/2), floor(ni/2))
    system.time(idx_p <- eval(parse(text=paste0("x1[-id,.I[order(-",var,")][1:",m_p,"]]"))))
    system.time(idx_n <- eval(parse(text=paste0("x1[-id,.I[order(",var,")][1:",ni-m_p,"]]"))))
    idx = c(idx_p,idx_n)
    id = c(id,idx)
    samp_info = rbind(samp_info,x1[idx,])
  }
  proc.time()-ptm
  return(samp_info)
}
