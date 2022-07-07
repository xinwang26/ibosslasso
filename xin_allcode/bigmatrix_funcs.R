#######################################################################
#Author: Xin Wang
#Initial Date: Tue May  3 20:33:58 2016
#Program Description: This script contains all functions needed for performing bid data computation which cannot be completed by existed functions in R base
#Completed:
#Dependencies:
#Main Reference:
#Input:
#Output: function sample_idx, getnorm,
#######################################################################

#######################################################################
#Function Description: this function returns the index to fit the regression
#Dependencies:
#Input: x -- all predictors
#       n_s -- number of observation to be kept
#       prop -- proportion of general norm account in the sample
#       sort_func -- the function used to perform the partial sorting
#Output: set of index of observation to be kept
#######################################################################
sample_idx <-function(x,n_s,prop = 0,sort_func = sort)
{
    n = nrow(x)
    p = ncol(x)
    k1 = floor(prop*n_s)#to include k largest norm observations
    id_taken = c()
    #first get the largest norm observations' indices
    if(prop>0){
        norm <- getnorm(x)
        id_taken <- which(norm>=(-sort(-norm, partial=k1)[k1]))
        remove(norm)
    }
    k2 = floor((n_s - k1)/p) #number of obs for each variable
    k21 = floor(k2/2)
    #print(k2)
    for(i in 1:p)
    {
        #print(length(id_taken))
        id_left = (1:n)[!(1:n)%in%id_taken]
        id_i= which(x[id_left,i]>=(-sort(-x[id_left,i],partial = k21)[k21]))#largest k2/2
        id_i1= which(x[id_left,i]<=(sort(x[id_left,i],partial = k21)[k21]))#smallest k2/2
        #sort works for big matrices
        id_taken = c(id_taken,id_i,id_i1)
    }
    return(id_taken)
}


#######################################################################
#Function Description:get norm for big matrix
#Dependencies: bigmemory, biganalytics packages
#Input:
#Output:
#######################################################################
getnorm<-function(bigmat,dim = 2){ #dim = 1 norm by row, dim = 2 norm by column
    #norm_entry <- apply(bigmat,dim,devi_sqr)#devi_sqr returns squared deviance from mean of corresponding dimension
    norm <- apply(as.big.matrix(apply(bigmat,dim,devi_sqr)),dim%%2+1,sum)
    gc()
    return(sqrt(norm))
}

#######################################################################
#Function Description: compute squared deviance from mean
#Dependencies: bigmemory, biganalytics
#Input:
#Output:
#######################################################################
devi_sqr <-function(bigvector){
    veclength = length(bigvector)
    return((bigvector - rep(mean(bigvector),veclength))^2)
}

#######################################################################
#Function Description: Lasso Sensitivity compute function
lasso_sens <- function(reg_v,true_v){sum(reg_v!=0&true_v!=0)/sum(reg_v!=0)}

#######################################################################
#Function Description: Lasso Specificity compute function
lasso_spec <- function(){sum(reg_v==0&true_v==0)/sum(reg_v==0)}

#######################################################################
#Function Description: Number of variable selected by Lasso
lasso_sele <- function(reg_v){sum(reg_v!=0)}
