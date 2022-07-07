
args <- commandArgs(trailingOnly = TRUE)
(n = as.numeric(args[1]))
(p = as.numeric(args[2]))
(subsize = as.numeric(args[3]))
r0 = subsize
test_size = 1000  #may need further adjustment
ifelse((exists("n") & !is.na(n)),print(n),n<-1e5)
ifelse((exists("p") & !is.na(p)),print(p),p<-1e3)
ifelse(exists("subsize")&(!is.na(subsize)),print(subsize),subsize<-1e4)
setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)
library(foreach)
library(doParallel)
registerDoParallel(10)  #parallel compute is faster but affect calculate of time cost
library(lars)#in case
library(glmnet)
source("data_Gen.R",echo = F)
#this contains functions for the ranking and selecting
source("data_table_funcs_fast.R",echo = F)
#library(irlba)
source("implement_module_fast_para.R",echo = F)
(datetime = substr(Sys.time(),6,19))

p_true = floor(sqrt(p)+1)
rounds = 1
set.seed(5)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,rounds)
seed2_list = sample(99999:199999,rounds)
seed3_list = sample(200000:299999,rounds)
seed4_list = sample(300000:399999,rounds)
seed4_list = sample(400000:499999,rounds)
seed5_list = sample(500000:599999,rounds)
seed6_list = sample(600000:699999,rounds)
seed7_list = sample(700000:799999,rounds)
(scale_beta  = sqrt(log(p)/subsize)*8)#necessary for data generation

X_distribution = "fast_t"
data_gen = data_gen_t
sigma1 = 1 #scale of error
source("evaluation.R",echo = F) # the evaluation function test-mse will rely on the distribution of predictors
#to store result
repo_s<- c()
ptm = proc.time()
convg_fail = data.frame()
cat("everything normal! before generating data")
source("implement_module_fast.R",echo=T)
for (i in 1:rounds){
  seed1 = seed1_list[i]#for x generation
  seed2 = seed2_list[i]#for beta generation
  seed3 = seed3_list[i]#for y generation
  seed4 = seed4_list[i]#for test x generation
  seed5 = seed5_list[i]#for test y generation
  seed6 = seed6_list[i]#for leverage sampling
  seed7 = seed7_list[i]#for logistic first sampling
  
  #generate the training data
  system.time(x <- data_gen(seed1))
  system.time(beta <- para_gen(seed2))
  true_beta = c(0,beta)
  system.time(x_mat <- as.matrix(x[,paste0("V",1:p),with=F]))
  z<-as.matrix(x_mat)%*%beta
  pr = 1/(1+exp(-z))
  set.seed(seed3)#for rbinom
  y <-  rbinom(dim(pr)[1],1,pr)  
  
  #generate the test data
  
  x_test <- data_gen(seed4,n_full = test_size)
  x_test <- as.matrix(x_test[,paste0("V",1:p),with=F])
  z_test <- as.matrix(x_test)%*%beta
  pr_test = 1/(1+exp(-z_test))
  set.seed(seed5)#for rbinom
  z_test <-  rbinom(dim(pr_test)[1],1,pr_test)  
  
  convg = 1
  #IBOSS subset and its lasso regression
  lasso_dopt = reg_module("inno",r=0)
  lasso_full = reg_module("full")
  convg = convg*(lasso_cor100$cvm_rch_min)
  
}


