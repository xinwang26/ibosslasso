#######################################################################
#Author: Xin Wang
#Initial Date: Tue Feb 14 01:13:36 2017
#Program Description: this program is designed to measure SIS's ability to include the true model by checking empirical distribution of s to include all true variables. Unlike other simulation, no need to consider compute time, only a limited number of dimension setting need to be considered
#Completed: No
#Dependencies: data_gen implement_model_fast  
#Main Reference: Sure independence screening by Fan
#Input: NA
#Output:
#######################################################################

args <- commandArgs(trailingOnly = TRUE)
(n = as.numeric(args[1]))
(p = as.numeric(args[2]))
(subsize = as.numeric(args[3]))

ifelse((exists("n") & !is.na(n)),print(n),n<-1e5)
ifelse((exists("p") & !is.na(p)),print(p),p<-1e3)
ifelse(exists("subsize")&(!is.na(subsize)),print(subsize),subsize<-1e4)
setwd("/home/wjin3/code/tocopy")
#install.packages("data.table")
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
source("implement_module_fast_para.R",echo = T)
(datetime = substr(Sys.time(),6,19))

p_true = floor(sqrt(p)+1)
rounds = 1000
set.seed(5)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,rounds)
seed2_list = sample(99999:199999,rounds)
seed3_list = sample(200000:299999,rounds)
seed4_list = sample(300000:399999,rounds)
seed4_list = sample(400000:499999,rounds)
seed5_list = sample(500000:599999,rounds)
seed6_list = sample(600000:699999,rounds)
(scale_beta  = sqrt(log(p)/subsize)/3)#necessary for data generation

X_distribution = "fast_t"
data_gen = data_gen_t
sigma1 = 1 #scale of error
source("evaluation.R",echo = F) # the evaluation function test-mse will rely on the distribution of predictors
#to store result
repo_s<- c()
ptm = proc.time()
convg_fail = data.frame()
cat("everything normal! before generating data")
conditions = data.frame("__",X_distribution,rounds,n,p,p_true,scale_beta,subsize)
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
  #x,x_mat,y,x_test,y_test will be used in all functions from now, passing by global environment
  # s_find<-function(pstar,p,...)
  
  tryCatch(s_value <- s_find(p_true,p,seed4),
           error = function(c){
             cat("the lasso/sortting cpp function has error, this",i,"-th iteration get skipped and s set as -1")
           s_value<- -1
           })
  repo_s <- c(repo_s,s_value)
  cat("the",i,"th run finished normally\n")
  write.table(repo_s,sep=",",file= paste0("s_find_result_from438to1000","n_",n,"ns_",subsize,"p_",p,"pstar_",p_true,".csv"))
  saveRDS(repo_s,file=paste0(substr(datetime,1,5),"_n_",n,"_p_",p,"_k_",subsize,"_d_",X_distribution,"sfind_",rounds,"_out.rds"))
  remove(x,y)
}

finalout = list(conditions,repo_s)
saveRDS(finalout,file=paste0(substr(datetime,1,5),"_n_",n,"_p_",p,"_k_",subsize,"_d_",X_distribution,"sfind_",rounds,"_finalout.rds"))




