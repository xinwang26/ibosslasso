setwd("/home/wjin3/code/tocopy")
library(data.table)
library(Metrics)

library(foreach)
library(doParallel)
registerDoParallel(10)
library(lars)
library(glmnet)
library(data.table)
library(glmnet)
#a = fread("./BlogFeedback/blogData_test-2012.02.01.00_00.csv",header=F)
source("data_table_funcs_fast.R",echo = F)
#source("implement_module_fast_para.R",echo = F)
# raw = fread("./BlogFeedback/blogData_train.csv") #import
raw = fread("YearPredictionMSD.txt")
colnames(raw) = colnames(raw)[c(91,1:90)]
head(raw)
raw  =cbind(raw[,c(2:91)],raw[,c(1)])
head(raw)

norm_data = cbind(scale(raw[,1:91])) #normalize
#norm_data[,c(13,33,38,278)] = raw[,c(13,33,38,278)]


#head(train_set)
result_compare = list()
n_selected = c()
out = list()
set.seed(2)
rounds=1000
seed_list = sample(1:999999,rounds*2)
#convg_fail = rbind(rep(0,(length(method_string)+1)),convg_fail)
MPEs = c()
MPE_calculator<- function(reg_list)
{
  (coef_l = coef(reg_list))
  y_hat = x_test%*%matrix(coef_l[-1],ncol=1)+ matrix(rep(1,dim(x_test)[1]),ncol =1)*coef_l[1]
  MPE = mse(y_hat,y_test)
  return(MPE)
}
timetrain = data.frame()
time_fulllist = c()
for(i in 1:rounds)
{
  set.seed(seed_list[i])
  train_index = c(1:463715)
  test_index = c(463716:nrow(raw))
  train_set = data.table(norm_data[train_index,])
  test_set = data.table(norm_data[test_index,])
  
  x_test = test_set[,c(1:90)]
  x_test = as.matrix(x_test)
  y_test = as.matrix(test_set[,91])
  
  
  names(train_set) = c(paste0("V",1: (dim(train_set)[2]-1)),"y")
  
  full_sample = train_set

  xtrain = full_sample[,c(1:90)]
  ytrain = as.matrix(full_sample[,c(91)])

  #coef(full_md)
  x_mat = as.matrix(xtrain)
  time_full = system.time(full_md <- cv.glmnet(x=x_mat,y=as.numeric(ytrain),alpha=1,family="gaussian"))
  time_fulllist = c(time_fulllist,time_full[3])
  time_dopt = system.time(dopt_subdata <- all_vars_max_real(x=xtrain,y = ytrain,sample_size=1800))
  time_cor = system.time(cor02_subdata <- cor_var_max_real(x=xtrain,y = ytrain,sample_size=1800,prop = 0.5))
  xtrain = xtrain[,c(1:90)]
  n_selected[i] = length(dopt_subdata)
  time_lev = system.time(lev_selec <- lev_sample_real(x = xtrain, y=ytrain,sample_size =n_selected[i],seed6=seed_list[i]))

  xtrain = as.matrix(xtrain)
  ytrain = as.matrix(full_sample[,91])
  time_dopt = time_dopt + system.time(dopt_train <- cv.glmnet(x=xtrain[c(dopt_subdata),],y=ytrain[c(dopt_subdata),],alpha=1,family="gaussian",parallel = T))
  time_cor = time_cor + system.time(cor02_train <- cv.glmnet(x=xtrain[c(cor02_subdata),],y=ytrain[c(cor02_subdata),],alpha=1,family="gaussian",parallel = ))
  time_lev = time_lev+ system.time(lev_train <- cv.glmnet(x=xtrain[c(lev_selec),],y=ytrain[c(lev_selec),],alpha=1,family="gaussian",parallel = T))
  set.seed(seed_list[rounds+1-i])
  time_rand = system.time(rand_selec <- sample(1:dim(full_sample)[1],n_selected[i]))
  time_rand = time_rand + system.time(rand_train <- cv.glmnet(x=xtrain[c(rand_selec),],y=ytrain[c(rand_selec),],alpha=1,family="gaussian",parallel = T))
  timetrain = rbind(timetrain,c(time_full,time_dopt,time_cor,time_lev,time_rand))
  result_compare[[i]] = cbind(FULL=coef(full_md),D_OPT=coef(dopt_train),COR02 = coef(cor02_train) ,LEV=coef(lev_train),UNIF=coef(rand_train))
  (nonzeros=unique(c(which(result_compare[[i]][,1]!=0),which(result_compare[[i]][,2]!=0),which(result_compare[[i]][,3]!=0),which(result_compare[[i]][,4]!=0),which(result_compare[[i]][,4]!=0) )))
  (out[[i]] = as.matrix(result_compare[[i]][sort(nonzeros),]))

      #raw_tests = fread(paste0("./BlogFeedback/blogData_test-2012.",month,".",day,".00_00.csv"))
  MPEs = rbind(MPEs,c(FULL = MPE_calculator(full_md),
                      COR02 = MPE_calculator(cor02_train),
                      D_OPT=MPE_calculator(dopt_train),
                      LEV=MPE_calculator(lev_train),
                      UNIF=MPE_calculator(rand_train)))

}
times
MPEs
saveRDS(out,"bootstrap_real1000_WITHCOR.rds")
saveRDS(result_compare,"bootstrap_real1000_full_WITHCOR.rds")

dist <-function(mat){
  (DOPT_dist = sum((mat[,1]-mat[,2])^2))
  (LEV_dist = sum((mat[,1]-mat[,3])^2))
  (UNIF_dist = sum((mat[,1]-mat[,4])^2))
  return(data.frame(DOPT_dist,LEV_dist,UNIF_dist))
}
aaa = lapply(out,dist)
bbb = unlist(aaa)
bbb = t(matrix(bbb,nrow=3))
apply(bbb,2,mean)
apply(bbb,2,sd)
ddd = apply(bbb,1,min)
sum(bbb[,1]==ddd)
sum(bbb[,2]==ddd)
sum(bbb[,3]==ddd)
ccc = unlist(unlist(result_compare))
saveRDS(result_compare,file="result_compare_1000.rds")
saveRDS(out,file="out_1000.rds")
