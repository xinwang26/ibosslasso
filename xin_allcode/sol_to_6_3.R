#######################################################################
#Author: Xin Wang
#Initial Date: Wed Jun 15 15:55:09 2016
#Program Description: when having 1e6*1e3 predictors, data is too large to use the original approach so here is solution for this specific case
#Completed: 
#Dependencies: 
#Main Reference: 
#Input: 
#Output:
#######################################################################

setwd("/Users/xinwang/Google Drive/research reading/2016Spring/BIGDATA")
library(bigmemory)
library(bigalgebra)
library(biglasso)
library(biganalytics)
library(data.table)
n = 1e6
p = 1e3
subsize = 1e4
n_true = floor(sqrt(p))
round =  1
seed1 = 100
seed2 = 555
seed3 = 999
scale_beta = sqrt(2*log(n)/subsize)
cat(sum(sapply(ls(),object.size)),"BYTES")
#system.time(x <- matrix(rnorm(n*p),n))
set.seed(seed1)
system.time(x <- as.data.table(matrix(rnorm(n*p),ncol=p)))
x[,id:=.I]
source("data_table_funcs.R")
ptm = proc.time()
informative_idx = all_vars_max(x,subsize) #this x need to have index in the first column
proc.time()-ptm
sub_id = as.numeric(informative_idx[,id])
sub_x <- as.big.matrix(informative_idx[,-1,with=F])
sub_y <- as.big.matrix(y[sub_id]) 
beep()
