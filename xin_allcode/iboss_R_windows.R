#######################################################################
#Author: Xin Wang
#Initial Date: Wed Apr 27 23:09:32 2016
#Program Description: this program perform ranking then regression with subset
#Completed: No
#Dependencies:
#Main Reference:
#Input:
#Output:
#######################################################################
setwd("C:/Users/Xin Wang/Google ‘∆∂À”≤≈Ã/research reading/2016Spring/BIGDATA")

library(bigmemory)
library(bigalgebra)
library(biglasso)
library(biganalytics)
library(data.table)
library(Metrics)

library(lars)
library(glmnet)
#this R script contains functions generating x and beta, argument only need two random seed and one value of scale for true parameter beta
source("data_Gen.R",echo = T,max.deparse.length = 10000)
#this contains functions for the ranking and selecting
source("data_table_funcs_small.R",echo = T,max.deparse.length = 10000)
#library(irlba)


datetime = substr(Sys.time(),6,16)
n = 1e4
p = 1e2
subsize = 1e3
test_size = 50
cor_prop = 0.05
n_true = floor(sqrt(p))
rounds = 100
set.seed(2)
seed1_list = sample(1:99999,rounds)
seed2_list = sample(99999:199999,rounds)
seed3_list = sample(200000:299999,rounds)
seed4_list = sample(300000:399999,rounds)
seed4_list = sample(400000:499999,rounds)
seed5_list = sample(500000:599999,rounds)
(scale_beta  = sqrt(log(p)/subsize/5))#necessary for data generation
sigma = 3



X_distribution = "t"
data_gen = data_gen_t
column_select = "corr"
row_select = cor_var_max #equally consider all variables / cor_var_max/ mod_var_max/cor_var_min/cor_to_max



source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors

sens = data.frame()
spec = data.frame()
mses = data.frame()
nums = data.frame()
sseb = data.frame()
ptm = proc.time()

for(i in 1:rounds)
{
  seed1 = seed1_list[i]
  seed2 = seed2_list[i]
  seed3 = seed3_list[i]
  seed4 = seed4_list[i]
  seed5 = seed5_list[i]

  system.time(x <- data_gen(seed1))#generate x and beta first 
  system.time(beta <- para_gen(seed2))
  true_beta = c(0,beta)


  system.time(x_mat <- as.big.matrix(x[,paste0("V",1:p),with=F]))
  set.seed(seed3)
  y <-  x_mat%*%beta+rnorm(n,0,sigma)
  
  system.time(informative_idx <- row_select(x,y,subsize,prop = cor_prop))
  sub_x <- as.big.matrix(x_mat[informative_idx,])
  sub_y = as.big.matrix(y[informative_idx])
  n_selected = length(informative_idx)

  set.seed(seed4)
  rand_idx = sample(1:dim(x)[1],n_selected) #this x need to have index in the first column
  proc.time()-ptm
  rand_x <- as.big.matrix(x_mat[rand_idx,])
  rand_y <-  as.big.matrix(y[rand_idx])

  system.time(sub_md_cv <- cv.biglasso(sub_x,sub_y,penalty = "lasso",family = "gaussian"))
  system.time(rand_md_cv <- cv.biglasso(rand_x,rand_y,penalty = "lasso",family = "gaussian"))
  system.time(md_cv <- cv.biglasso(x_mat,y,penalty = "lasso",family = "gaussian"))

  nums = rbind(nums,c(num_catch(coef(md_cv)),
                      num_catch(coef(sub_md_cv)),
                      num_catch(coef(rand_md_cv))))
  sens = rbind(sens,c(sensitivity(coef(md_cv)),
                      sensitivity(coef(sub_md_cv)),
                      sensitivity(coef(rand_md_cv))))
  spec = rbind(spec,c(specificity(coef(md_cv)),
                      specificity(coef(sub_md_cv)),
                      specificity(coef(rand_md_cv))))
  temses = rbind(mses,c(test_mse(md_cv,test_size,seed5,sigma1 = sigma),
                      test_mse(sub_md_cv,test_size,seed5,sigma1 = sigma),
                      test_mse(rand_md_cv,test_size,seed5,sigma1 = sigma)))
  trmses = rbind(mses,c(train_mse(x_mat,y,md_cv),
                        train_mse(x_mat,y,sub_md_cv),
                        train_mse(x_mat,y,rand_md_cv)))
  sseb = rbind(sseb,c(sse_beta(coef(md_cv)),
                       sse_beta(coef(sub_md_cv)),
                       sse_beta(coef(rand_md_cv))))

}

(output_names = c("full",paste(X_distribution,column_select,"informative",sep = "-"),"random"))

names(spec) = output_names
names(sens) = output_names
names(trmses) = output_names
names(temses) = output_names
names(sseb) = output_names

#write.csv(sens,file = paste0("sensi",datetime,".csv"))
#write.csv(spec,file = paste0("spec",datetime,".csv"))
#write.csv(temses,file = paste0("trmses",datetime,".csv"))
#write.csv(trmses,file = paste0("trmses",datetime,".csv"))
#write.csv(sseb,file = paste0("mses",datetime,".csv"))
sink(file=paste(substr(datetime,1,8),"result.txt"))
cat("data size:")
n
cat("number of variable:")
p
cat("number of most correlated variables optimized:")
cor_prop
cat("error term standard deviation is:")
sigma
cat("intended sample size:")
subsize
cat("implemented sample size:")
n_selected
cat("number of true variable:")
n_true
cat("simulation times:")
rounds
cat("scale of true parameters:")
scale_beta#necessary for data genera
cat("\nsensitivity:\n")
apply(sens,2,mean)#larger better
cat("specificity:\n")
apply(spec,2,mean)#larger better
cat("training mses:\n")
apply(trmses,2,mean)#smaller better
cat("testing mses:\n")
apply(temses,2,mean)#smaller better
cat("number of variable catched:\n")
apply(nums,2,mean)#smaller better
cat("total squared estimator error from true-values:\n")
apply(sseb,2,mean)#smaller better
cat("running time:\n")
proc.time()-ptm
sink()




