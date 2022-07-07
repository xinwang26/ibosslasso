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
#n = 10000
p = 5000
subsize = 5000
#args <- commandArgs(trailingOnly = TRUE)
#(n = as.numeric(args[1]))
#(p = as.numeric(args[2]))
#(subsize = as.numeric(args[3]))
setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)
library(foreach)
library(doParallel)
registerDoParallel(10)
library(lars)
library(glmnet)
source("data_Gen.R",echo = F)
#this contains functions for the ranking and selecting
source("data_table_funcs_fast.R",echo = F)
#library(irlba)
source("implement_module_fast_para.R",echo = F)
#library(irlba)
(datetime = substr(Sys.time(),6,19))

test_size = 1000
cor_prop = 0.2
p_true = floor(sqrt(p)+1)
rounds = 100
set.seed(5)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,rounds)
seed2_list = sample(99999:199999,rounds)
seed3_list = sample(200000:299999,rounds)
seed4_list = sample(300000:399999,rounds)
seed4_list = sample(400000:499999,rounds)
seed5_list = sample(500000:599999,rounds)
seed6_list = sample(600000:699999,rounds)
(scale_beta  = sqrt(log(p)/1000)/1.5)#necessary for data generation

X_distribution = "fast_t"
data_gen = data_gen_t
sigma1 = 1 #scale of error
source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors
#to store result
#repo_all = data.frame()
repo_cor100 = data.frame()
#repo_cor10 = data.frame()
repo_cor50 = data.frame()
repo_rand = data.frame()
#repo_lev = data.frame()
#repo_full = data.frame()
ptm = proc.time()
convg_fail = data.frame()

cat("everything normal!")
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
  convg = 1
  #IBOSS subset and its lasso regression
  lasso_cor100 = reg_module("cor",r=100/p)
  convg = convg*(lasso_cor100$cvm_rch_min)
  #lasso_cor10 = reg_module("cor",r=10/p)
  lasso_cor50 = reg_module("cor",r=50/p)
  convg = convg*(lasso_cor50$cvm_rch_min)
  #lasso_full = reg_module("full")
  #convg = convg*(lasso_full$cvm_rch_min)
  lasso_random = reg_module("rand")
  convg = convg*(lasso_random$cvm_rch_min)
  
  if(convg==1){
    
    #row_all = eval_module(lasso_all)
    row_cor100 = eval_module(lasso_cor100)
    #row_cor10 = eval_module(lasso_cor10)
    row_cor50 = eval_module(lasso_cor50)
    #row_lev = eval_module(lasso_lev)
    #row_full = eval_module(lasso_full)
    row_rand = eval_module(lasso_random)
    
    #repo_all = rbind(repo_all,row_all)
    repo_cor100 = rbind(repo_cor100,row_cor100)
    #repo_cor10 = rbind(repo_cor10,row_cor10)
    repo_cor50 = rbind(repo_cor50,row_cor50)
    #repo_lev = rbind(repo_lev,row_lev)
    #repo_full = rbind(repo_full,row_full)
    repo_rand = rbind(repo_rand,row_rand)
  }
  if(convg==0){
    rowi = c(i,
             #lasso_full$cvm_rch_min,
             lasso_cor100$cvm_rch_min,
             lasso_cor50$cvm_rch_min,
             #lasso_cor10$cvm_rch_min,
             lasso_random$cvm_rch_min)
    convg_fail = rbind(convg_fail,rowi)
  }
  
  
}

method_string = c(#"full",
                  "cor100",
                  "cor50",
                 # "cor10",
                  #"lev",
                  "rand")
names(convg_fail) = c("no_run",method_string)
output = data.frame()
for (m in method_string)
{
  output = rbind(output,
                 reportmean_module(get(paste0("repo_",m)))
  )}
rownames(output) = method_string
total_time = (proc.time()- ptm)[3]
conditions = data.frame("__",X_distribution,rounds,n,p,p_true,scale_beta,subsize,total_time,test_size)
print(conditions)
print(output)
cat("total time is:",total_time)
#output to csv
write.table("",sep=",",file= paste0("reg_R_result_table",datetime,".csv"))
write.table(conditions,sep=",",file = paste0("reg_R_result_table",datetime,".csv"),append=T,col.names = T,row.names = T)
write.table(output,sep = ",",file = paste0("reg_R_result_table",datetime,".csv"),append = T,col.names = T,row.names = T)
for (m in method_string)
{
  write.table(conditions,sep=",",file = paste0(m,"_",gsub("[[:blank:]]", "", datetime),".csv"),append=T,col.names = T,row.names = T)
  write.table(get(paste0("repo_",m)),sep = ",",file = paste0(m,"_",gsub("[[:blank:]]", "", datetime),".csv"),append = T,col.names = T,row.names = T)
  }

write.table(conditions,sep = ",",file = paste0("converge_",datetime,".csv"),append = T,col.names = T,row.names = T)
write.table(convg_fail,sep = ",",file = paste0("converge_",datetime,".csv"),append = T,col.names = T,row.names = T)