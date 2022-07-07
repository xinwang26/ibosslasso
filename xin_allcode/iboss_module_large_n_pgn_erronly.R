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
ifelse((exists("p") & !is.na(p)),print(p),p<-5000)
ifelse(exists("subsize")&(!is.na(subsize)),print(subsize),subsize<-2e3)
setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)
library(foreach)
library(doParallel)
registerDoParallel(1)  #parallel compute is faster but affect calculate of time cost
if(n*p>1e9){registerDoParallel(1) }
library(lars)
library(glmnet)
source("data_Gen.R",echo = F)
#this contains functions for the ranking and selecting
source("data_table_funcs_fast.R",echo = F)
#library(irlba)
source("implement_module_fast_nonpa.R",echo = F)
#library(irlba)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<3){print("wrong input!")} 
(datetime = substr(Sys.time(),6,19))

test_size = 1000
cor_prop = 0.2
p_true = floor(sqrt(p)+1)
rounds = 100
set.seed(3)
islasso = 1 #parameter for lasso or ridge or glmnet
seed1_list = sample(1:99999,5000)
seed2_list = sample(99999:199999,5000)
seed3_list = sample(200000:299999,5000)
seed4_list = sample(300000:399999,5000)
seed4_list = sample(400000:499999,5000)
seed5_list = sample(500000:599999,5000)
seed6_list = sample(600000:699999,5000)
method_string = c("full"
                  ,"cor01"
                  #            ,"cor02"
                  ,"cor001"
                  ,"cor002"
                  ,"cor005"
                  ,"lev"
                  ,"rand"
)
(scale_beta  = sqrt(log(p)/subsize)*1.5)#necessary for data generation

X_distribution = "fast_n"
data_gen = data_gen_n
sigma1 = 1 #scale of error
source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors
#to store result
repo_all = data.frame()
repo_cor01 = data.frame()
repo_cor02 = data.frame()
repo_cor001 = data.frame()
repo_cor002 = data.frame()
repo_cor005 = data.frame()
repo_rand = data.frame()
repo_lev = data.frame()
repo_full = data.frame()
convg_fail = data.frame()
convg_fail = rbind(convg_fail,c("nbr",rep(0,9)))
ptm = proc.time()
conditions = data.frame("__",X_distribution,rounds,n,p,p_true,scale_beta,subsize,test_size)
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
  lasso_cor01 = reg_module("cor",r=floor(0.1*p)/p)
  convg = convg*(lasso_cor01$cvm_rch_min)
  
  #lasso_cor02 = reg_module("cor",r=floor(0.2*p)/p)
  #convg = convg*(lasso_cor02$cvm_rch_min)
  
  lasso_cor001 = reg_module("cor",r=floor(0.01*p)/p)
  convg = convg*(lasso_cor001$cvm_rch_min)
  
  lasso_cor002 = reg_module("cor",r=floor(0.02*p)/p)
  convg = convg*(lasso_cor001$cvm_rch_min)
  
  lasso_cor005 = reg_module("cor",r=floor(0.05*p)/p)
  convg = convg*(lasso_cor001$cvm_rch_min)
  
  lasso_full = reg_module("full",r=0)
  convg = convg*(lasso_full$cvm_rch_min)
  
  lasso_lev = reg_module("lev",r=0)
  convg = convg*(lasso_lev$cvm_rch_min)
  
  lasso_random = reg_module("rand",r=0)
  convg = convg*(lasso_random$cvm_rch_min)
  
  if(convg==1){
    row_cor01 = eval_module(lasso_cor01)
   # row_cor02 = eval_module(lasso_cor02)
    row_cor001 = eval_module(lasso_cor001)
    row_cor002 = eval_module(lasso_cor002)
    row_cor005 = eval_module(lasso_cor005)
    row_lev = eval_module(lasso_lev)
    row_full = eval_module(lasso_full)
    row_rand = eval_module(lasso_random)
    
    repo_cor01 = rbind(repo_cor01,row_cor01)
    #repo_cor02 = rbind(repo_cor02,row_cor02)
    repo_cor001 = rbind(repo_cor001,row_cor001)
    repo_cor002 = rbind(repo_cor002,row_cor002)
    repo_cor005 = rbind(repo_cor005,row_cor005)
    repo_lev = rbind(repo_lev,row_lev)
    repo_full = rbind(repo_full,row_full)
    repo_rand = rbind(repo_rand,row_rand)
  }
  if(convg==0){
    rounds = rounds +1
    rowi = c(i,
             lasso_full$cvm_rch_min,
             lasso_cor01$cvm_rch_min,
     #        lasso_cor02$cvm_rch_min,
             lasso_cor001$cvm_rch_min,
             lasso_cor002$cvm_rch_min,
             lasso_cor005$cvm_rch_min,
             lasso_lev$cvm_rch_min,
             lasso_random$cvm_rch_min
   )
    convg_fail = rbind(convg_fail,rowi)
  }
}

output = data.frame()
names(convg_fail) = c("nBR_run",method_string)
failure_times = dim(convg_fail)[1]-1
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
