#######################################################################
#Author: Xin Wang
#Initial Date: Mon Sep  4 15:25:28 2017
#Program Description: 
#Completed: 
#Dependencies: this program implement the split and conquer approach 
#Main Reference: 
#Input: 
#Output:
#######################################################################
rounds = 100
args <- commandArgs(trailingOnly = TRUE)
(n = as.numeric(args[1]))
(p = as.numeric(args[2]))
(subsize = as.numeric(args[3]))
ifelse((exists("n") & !is.na(n)),print(n),n<-1e4)
ifelse((exists("p") & !is.na(p)),print(p),p<-5e3)
ifelse(exists("subsize")&(!is.na(subsize)),print(subsize),subsize<-1e3)
setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)
library(foreach)
library(doParallel)
library(lars)
library(glmnet)
source("data_Gen.R",echo = F)
#this contains functions for the ranking and selecting
source("data_table_funcs_fast.R",echo = F)
#library(irlba)
source("implement_module_fast_para.R",echo = F)

#library(irlba)
if(length(args)<3){print("wrong input!")} 
(datetime = substr(Sys.time(),6,19))
lam_min_rate = 0.005
num_slice = 10
omega = 2
test_size = 1000
cor_prop = 0.1
if(p==5000){(p_true = floor(sqrt(p)+1))}
if(p!=5000){(p_true = round(73/5000*p,0))}

set.seed(5)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,5000)
seed2_list = sample(99999:199999,5000)
seed3_list = sample(200000:299999,5000)
seed4_list = sample(300000:399999,5000)
seed4_list = sample(400000:499999,5000)
seed5_list = sample(500000:599999,5000)
seed6_list = sample(600000:699999,5000)
(scale_beta  = sqrt(log(p)/1000)/2)#necessary for data generation

X_distribution = "fast_t"
data_gen = data_gen_t
sigma1 = 1 #scale of error
source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors
#to store result
# repo_all = data.frame()
#repo_cor002 = data.frame()
#repo_cor002 = data.frame()
#repo_cor001 = data.frame()
repo_corsis005 = data.frame()
repo_randsis005 = data.frame()
#repo_lev = data.frame()
repo_sis005spc2o5 = data.frame()
repo_sis005spc4o10 = data.frame()
repo_sis005spc45o20 = data.frame()
repo_sisfull = data.frame()
ptm = proc.time()
convg_fail = data.frame()
method_string = c(
  #"full",
  # "all",
  #"lev",
  #"cor001",
  "sisfull",
  'sis005spc2o5',
  "sis005spc4o10",
  "sis005spc45o20",
  "corsis005",
  #"cor002",
  #"lev",
  "randsis005")
convg_fail = rbind(rep(0,(length(method_string)+1)),convg_fail)
registerDoParallel(1)
cat("everything normal!")
for(i in 1:rounds)
{
  
  #Data generate (makesure same random seed kept)
  #######################################################################
  cat("\nthe",i,"th round start")
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

  convg = 1
  lasso_sisfull = sis("full",cor_r=0.05)
  convg = convg*(lasso_sisfull$cvm_rch_min)
  
  lasso_sis005spc2o5 = sis("spc",cor_r=0.4,sis_r=0.05,num_slice = 5)
  convg = convg*(lasso_sis005spc2o5$cvm_rch_min)
  
  lasso_sis005spc4o10 = sis("spc",cor_r=0.4,sis_r=0.05,num_slice = 10)
  convg = convg*(lasso_sis005spc4o10$cvm_rch_min)
  
  lasso_sis005spc45o20 = sis("spc",cor_r=0.45,sis_r=0.05,num_slice = 20)
  convg = convg*(lasso_sis005spc45o20$cvm_rch_min)  

  lasso_corsis005 = sis("cor",cor_r=0.05)
  convg = convg*(lasso_corsis005$cvm_rch_min)
  
  #lasso_cor001 = reg_module("cor",cor_r=0.01)
  #convg = convg*(lasso_cor001$cvm_rch_min)
  
  #lasso_cor002 = reg_module("cor",cor_r=0.02)
  #convg = convg*(lasso_cor002$cvm_rch_min)
  
  lasso_randsis005 = sis("rand",cor_r=0.05)
  convg = convg*(lasso_randsis005$cvm_rch_min)

  row_corsis005 = eval_module(lasso_corsis005)
  # row_cor001 = eval_module(lasso_cor001)
  # row_cor002 = eval_module(lasso_cor002)
  row_sisfull= eval_module(lasso_sisfull)
  row_sis005spc2o5 = eval_module(lasso_sis005spc2o5)
  row_sis005spc4o10 = eval_module(lasso_sis005spc4o10)
  row_sis005spc45o20 = eval_module(lasso_sis005spc45o20)
  row_randsis005 = eval_module(lasso_randsis005)
  
  repo_corsis005 = rbind(repo_corsis005,row_corsis005)
  # repo_cor001 = rbind(repo_cor001,row_cor001)
  # repo_cor002 = rbind(repo_cor002,row_cor002)
  
  repo_sis005spc2o5 = rbind(repo_sis005spc2o5,row_sis005spc2o5)
  repo_sis005spc4o10 = rbind(repo_sis005spc4o10,row_sis005spc4o10)
  repo_sis005spc45o20 = rbind(repo_sis005spc45o20,row_sis005spc45o20)
  repo_randsis005 = rbind(repo_randsis005,row_randsis005)
  repo_sisfull = rbind(repo_sisfull,row_sisfull)
  if(convg==0){
    #rounds = rounds +1
    rowi = c(i,
             lasso_sisfull$cvm_rch_min,
             lasso_sis005spc2o5$cvm_rch_min,
             lasso_sis005spc4o10$cvm_rch_min,
             lasso_sis005spc45o20$cvm_rch_min,
             lasso_corsis005$cvm_rch_min,
             # lasso_cor001$cvm_rch_min,
             # lasso_cor002$cvm_rch_min,
             lasso_randsis005$cvm_rch_min)
    convg_fail = rbind(convg_fail,rowi)
  }
}


names(convg_fail) = c("no_run",method_string)
output = data.frame()
rd_repo = list()
rd_record = list()
for (m in method_string)
{
  rd_repo[[m]] = round(apply(get(paste0("repo_",m)),2,mean),4)
  rd_record[[m]] = get(paste0("repo_",m))
}

total_time = (proc.time()- ptm)[3]
conditions = data.frame("__",X_distribution,rounds,n,p,p_true,scale_beta,subsize,total_time,test_size)
rd_repo[['conditions']] = conditions
rd_repo[['converge']] = convg_fail

print(rd_repo)


saveRDS(rd_repo,file=paste0(substr(datetime,1,5),"_n_",n,"_p_",p,"_k_",subsize,"_d_",X_distribution,"runs_",rounds,"_splitconq","_out.rds"))
saveRDS(rd_record,file=paste0(substr(datetime,1,5),"_n_",n,"_p_",p,"_k_",subsize,"_d_",X_distribution,"runs_",rounds,"_splitconq","_record.rds"))



