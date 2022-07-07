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
args <- commandArgs(trailingOnly = TRUE)
(n = as.numeric(args[1]))
(p = as.numeric(args[2]))
(subsize = as.numeric(args[3]))
ifelse((exists("n") & !is.na(n)),print(n),n<-1e5)
ifelse((exists("p") & !is.na(p)),print(p),p<-5e3)
ifelse(exists("subsize")&(!is.na(subsize)),print(subsize),subsize<-3e3)
setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)
library(foreach)
library(doParallel)
registerDoParallel(5)
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

test_size = 1000
cor_prop = 0.1
p_true = floor(sqrt(p)+1)
rounds = 100
set.seed(5)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,5000)
seed2_list = sample(99999:199999,5000)
seed3_list = sample(200000:299999,5000)
seed4_list = sample(300000:399999,5000)
seed4_list = sample(400000:499999,5000)
seed5_list = sample(500000:599999,5000)
seed6_list = sample(600000:699999,5000)
(scale_beta  = sqrt(log(p)/1000))#necessary for data generation

X_distribution = "fast_m"
data_gen = data_gen_m
sigma1 = 1 #scale of error
source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors
#to store result
#repo_all = data.frame()
repo_cor002 = data.frame()
repo_cor01 = data.frame()
repo_cor001 = data.frame()
repo_rand = data.frame()
#repo_lev = data.frame()
repo_full = data.frame()
ptm = proc.time()
convg_fail = data.frame()
method_string = c(
  "full",
 # "all",
 # "lev",
  "cor002",
  "cor001",
  "cor01",
  #"lev",
  "rand")
convg_fail = rbind(rep(0,(length(method_string)+1)),convg_fail)

cat("everything normal!")
for(i in 1:rounds)
{
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
  #x,x_mat,y,x_test,y_test will be used in all functions from now, passing by global environment
  
  #IBOSS subset and its lasso regression
  #lasso_all = reg_module("all")
  convg = 1
  #IBOSS subset and its lasso regression
  lasso_cor002 = reg_module("cor",r=0.02)
  convg = convg*(lasso_cor002$cvm_rch_min)
  
  lasso_cor01 = reg_module("cor",r=0.1)
  convg = convg*(lasso_cor01$cvm_rch_min)
  
  lasso_cor001 = reg_module("cor",r=0.01)
  convg = convg*(lasso_cor001$cvm_rch_min)
  
 # lasso_all = reg_module("all")
#  convg = convg*(lasso_all$cvm_rch_min)
  
 # lasso_lev = reg_module("lev")
#  convg = convg*(lasso_lev$cvm_rch_min)
  
  lasso_full = reg_module("full")
  convg = convg*(lasso_full$cvm_rch_min)
  
  lasso_random = reg_module("rand")
  convg = convg*(lasso_random$cvm_rch_min)
  
  
#  row_all = eval_module(lasso_all)
  row_cor002 = eval_module(lasso_cor002)
  row_cor01 = eval_module(lasso_cor01)
  row_cor001 = eval_module(lasso_cor001)
  #row_lev = eval_module(lasso_lev)
  row_full = eval_module(lasso_full)
  row_rand = eval_module(lasso_random)
  
  #repo_all = rbind(repo_all,row_all)
  repo_cor002 = rbind(repo_cor002,row_cor002)
  repo_cor01 = rbind(repo_cor01,row_cor01)
  repo_cor001 = rbind(repo_cor001,row_cor001)
 # repo_lev = rbind(repo_lev,row_lev)
  repo_full = rbind(repo_full,row_full)
  repo_rand = rbind(repo_rand,row_rand)
  
  if(convg==0){
    rounds = rounds +1
    rowi = c(i,
             lasso_full$cvm_rch_min,
             #lasso_all$cvm_rch_min,
            # lasso_lev$cvm_rch_min,
             lasso_cor002$cvm_rch_min,
             lasso_cor001$cvm_rch_min,
             lasso_cor01$cvm_rch_min,
             lasso_random$cvm_rch_min)
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

saveRDS(rd_repo
        ,file=paste0(substr(datetime,1,5),"_n_",n,"_p_",p,"_k_",subsize,"_d_",X_distribution,"_out.rds"))
saveRDS(rd_record
        ,file=paste0(substr(datetime,1,5),"_n_",n,"_p_",p,"_k_",subsize,"_d_",X_distribution,"_record.rds"))


print(rd_repo)
