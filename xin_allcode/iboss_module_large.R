#######################################################################
#Author: Xin Wang
#Initial Date: Wed Apr 27 23:09:32 2016
#Program Description: this program perform ranking then regression with subset, using all regular R data structure and regular lasso package (not big matrix)
#Completed: No
#Dependencies:
#Main Reference:
#Input:
#Output:
#######################################################################

setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)
library(foreach)
library(doParallel)
registerDoParallel(4)
library(lars)
library(glmnet)
source("data_Gen.R",echo = T,max.deparse.length = 10000)
#this contains functions for the ranking and selecting
source("data_table_funcs_fast.R",echo = T,max.deparse.length = 10000)
#library(irlba)
source("implement_module_fast.R",echo = T,max.deparse.length = 10000)
#library(irlba)

datetime = substr(Sys.time(),6,16)
n = 1e5
p = 1e3
subsize = 1e4
test_size = n/100
cor_prop = 0.2
n_true = floor(sqrt(p)+1)
rounds = 100
set.seed(3)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,rounds)
seed2_list = sample(99999:199999,rounds)
seed3_list = sample(200000:299999,rounds)
seed4_list = sample(300000:399999,rounds)
seed4_list = sample(400000:499999,rounds)
seed5_list = sample(500000:599999,rounds)
seed6_list = sample(600000:699999,rounds)
(scale_beta  = 1.5*sqrt(log(p)/subsize))#necessary for data generation

X_distribution = "fast_n"
data_gen = data_gen_n
sigma1 = 1 #scale of error
source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors
#to store result
repo_all = data.frame()
repo_cor01 = data.frame()
repo_cor02 = data.frame()
repo_cor03 = data.frame()
repo_rand = data.frame()
repo_lev = data.frame()
repo_full = data.frame()
ptm = proc.time()

for(i in 1:rounds)
{
  seed1 = seed1_list[i]#for x generation
  seed2 = seed2_list[i]#for beta generation
  seed3 = seed3_list[i]#for y generation
  seed4 = seed4_list[i]#for test x generation
  seed5 = seed5_list[i]#for test y generation
  seed6 = seed6_list[i]#for leverage sampling
  
  #generate the training data
  system.time(x <- data_gen(seed1))
  system.time(beta <- para_gen(seed2))
  true_beta = c(0,beta)
  system.time(x_mat <- as.matrix(x[,paste0("V",1:p),with=F]))
  set.seed(seed3)
  y <-  as.matrix(x_mat)%*%beta+rnorm(n,0,sigma1)
  
  #generate the test data
  x_test <- data_gen(seed4,n_full = test_size)
  x_test <- as.matrix(x_test[,paste0("V",1:p),with=F])
  set.seed(seed5)
  y_test <- x_test%*%beta+rnorm(test_size,0,sigma1) 
  y_true <- x_test%*%beta  
  #x,x_mat,y,x_test,y_test will be used in all functions from now, passing by global environment
  
  #IBOSS subset and its lasso regression
  lasso_all = reg_module("all")
  lasso_cor02 = reg_module("cor",r=0.2)
  lasso_cor03 = reg_module("cor",r=0.3)
  lasso_cor01 = reg_module("cor",r=0.1)
  lasso_lev = reg_module("lev",seed6)
  lasso_full = reg_module("full")
  lasso_random = reg_module("rand")
  
  
  row_all = eval_module(lasso_all)
  row_cor01 = eval_module(lasso_cor01)
  row_cor02 = eval_module(lasso_cor02)
  row_cor03 = eval_module(lasso_cor03)
  row_lev = eval_module(lasso_lev)
  row_full = eval_module(lasso_full)
  row_rand = eval_module(lasso_random)
  
  repo_all = rbind(repo_all,row_all)
  repo_cor01 = rbind(repo_cor01,row_cor01)
  repo_cor02 = rbind(repo_cor02,row_cor02)
  repo_cor03 = rbind(repo_cor03,row_cor03)
  repo_lev = rbind(repo_lev,row_lev)
  repo_full = rbind(repo_full,row_full)
  repo_rand = rbind(repo_rand,row_rand)
  
}

method_string = c("full","cor01","cor02","cor03","all","lev","rand")
output = data.frame()
for (m in method_string)
{
  output = rbind(output,
                 reportmean_module(get(paste0("repo_",m)))
  )}
rownames(output) = method_string
total_time = (proc.time()- ptm)[3]
conditions = data.frame("__",X_distribution,rounds,n,p,n_true,scale_beta,subsize,total_time)
print(conditions)
print(output)
cat("total time is:",total_time)
#output to csv
write.table("",sep=",",file= paste0("reg_R_result_table",datetime,".csv"))
write.table(conditions,sep=",",file = paste0("reg_R_result_table",datetime,".csv"),append=T,col.names = T,row.names = T)
write.table(output,sep = ",",file = paste0("reg_R_result_table",datetime,".csv"),append = T,col.names = T,row.names = T)

