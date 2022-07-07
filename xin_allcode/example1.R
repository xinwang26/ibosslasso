setwd("/home/wjin3/code/tocopy")
library(data.table)
library(glmnet)
#a = fread("./BlogFeedback/blogData_test-2012.02.01.00_00.csv",header=F)
source("data_table_funcs_fast.R",echo = F)
#source("implement_module_fast_para.R",echo = F)
raw = fread("./BlogFeedback/blogData_train.csv") #import
test1 = fread("./BlogFeedback/blogData_test-2012.03.17.00_00.csv") #import

#train_set = data.table(train_set)
#train_set[,c(13,33,38,278)] = raw[,c(13,33,38,278)]

pre_proc <- function(data){
  data1 = cbind(scale(data[,1:280]),data[,281])
  data1 = data.table(data1)
  #data1[,c(13,33,38,278)] = data[,c(13,33,38,278)]
  nacols = c()
  for(coli in 1:ncol(data)){
    if(is.na(data1[1,coli,with=F])){
      nacols = c(nacols,coli)
    }}
  
  nacols
  data1[,c(nacols)] = data[,c(nacols),with=F]
  return(data1)
}
train_set = pre_proc(raw) #normalize
test1 = pre_proc(test1)

names(train_set) = c(paste0("V",1: (dim(train_set)[2]-1)),"y")
#head(train_set)
result_compare = list()
n_selected = c()
out = list()
set.seed(2)
rounds=1000
seed_list = sample(1:999999,rounds*2)
i=1
predmse<-function(data,model){
  names(data)[ncol(data)] ="y" 
  y = data[,"y"]
  x = data[,-ncol(data),with=F]
  predictions = predict(model,as.matrix(x))
  mse_val = sum((y-predictions)^2)/nrow(data)
  return(mse_val)
  }
for(i in 1:rounds)
{
  set.seed(seed_list[iter])
  full_sample = train_set
  #[sample(1:dim(train_set)[1],5e4,replace = T),]
  
  xtrain = full_sample[,c(1:280)]
  ytrain = as.matrix(full_sample[,c(281)])
  
  #coef(full_md)
  x_mat = as.matrix(xtrain)
  full_md = cv.glmnet(x=x_mat,y=as.numeric(ytrain),alpha=1,family="gaussian")
  dopt_subdata =all_vars_max_real(x=xtrain,y = ytrain,sample_size=5600)
  dopt_subdata =cor_var_max_real(x=xtrain,y = ytrain,sample_size=5600,prop=0.2)
  xtrain = xtrain[,c(1:280)]
  n_selected[iter] = length(dopt_subdata)
  lev_selec = lev_sample_real(x = xtrain, y=ytrain,sample_size =n_selected[iter],seed6=seed_list[iter])
  
  xtrain = as.matrix(xtrain)
  ytrain = as.matrix(full_sample[,281])
  dopt_train = cv.glmnet(x=xtrain[c(dopt_subdata),],y=ytrain[c(dopt_subdata),],alpha=1,family="gaussian")
  lev_train = cv.glmnet(x=xtrain[c(lev_selec),],y=ytrain[c(lev_selec),],alpha=1,family="gaussian")
  set.seed(seed_list[iter])
  rand_selec = sample(1:dim(full_sample)[1],n_selected[iter])
  rand_train = cv.glmnet(x=xtrain[c(rand_selec),],y=ytrain[c(rand_selec),],alpha=1,family="gaussian")
  
  
  #################################
  #split and conquer
  num_slice = 10
  vote_rate =0.45
  omega = floor(vote_rate*num_slice)
  n=50000
  nk = floor(n/num_slice)  #num_slice is just k in the paper, name k is probably occupied, incase override use num_slice
  n_selected = nk
  nk_end = n - (num_slice-1)*nk #make sure the last share cover the remain
  y=ytrain
  library(doParallel)
  registerDoParallel(num_slice)
  dist_resp= "gaussian"
  islasso=1
  split_md <- foreach( i = 1:num_slice
                       
                       # ,.multicombine=TRUE
                       # ,.init=list(list(), list()) #combine result as list will make it easier to combine?
  ) %dopar% 
  {
    n_k= ifelse(i<num_slice,nk,nk_end)
    slice_start = (i-1)*nk+1
    slice_end = (i-1)*nk+n_k
    #fit lasso
    slice_x <- x_mat[slice_start:slice_end,]
    slice_y <-  y[slice_start:slice_end]
    md_cv <- cv.glmnet(slice_x,slice_y,parallel = T,alpha=islasso,family = dist_resp)
    cvm_rch_min = (md_cv$lambda.min!=md_cv$lambda[length(md_cv$lambda)])
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    betak=list(betahat=coef(md_cv),mse.min=mse.min,cvm_rch_min = 1*cvm_rch_min)
  }
  #lasso_time = proc.time()-ptm_start
  coef_list = lapply(split_md, `[[`, 1)
  split_beta = do.call(cbind, unlist(coef_list, recursive=TRUE))
  mse.min =  mean(unlist(lapply(split_md, `[[`, 2), recursive=TRUE))
  cvm_rch_min = min(unlist(lapply(split_md, `[[`, 3), recursive=TRUE))
  #combine
  #######################################################################
  theta = list()
  nonzero_ind= as.matrix(1*(split_beta!=0))
  #tic = proc.time()
  Ac = apply(nonzero_ind,1,sum)
  Ac = 1* (Ac > omega)
  indicates = which(Ac!=0)
  Ac1 = Ac[indicates]
  E =  as(diag(Ac1), "dgCMatrix") #this is A
  component = list()
  betai = list()
  Sigma = list()
  for (i in 1:num_slice){
    n_k= ifelse(i<num_slice,nk,nk_end)
    slice_start = (i-1)*nk+1
    slice_end = (i-1)*nk+n_k
    Xk = as.matrix(cbind(rep(1,n_k), x_mat[slice_start:slice_end,] )[,indicates])
    system.time(theta[[i]] <- Xk %*% split_beta[indicates,i] )
    system.time(component[[i]] <- E%*%t(Xk)%*%Xk%*%E)
    system.time(betai[[i]] <- component[[i]]%*%split_beta[indicates,i])
  }
  system.time(beta_combine <- E %*% ginv(as.matrix(Reduce(`+`, component))) %*% Reduce(`+`, betai))
  #comb_time = proc.time()-tic
  p=280
  beta_restore = as(matrix(rep(0,p+1),ncol=1),"dgCMatrix")
  beta_restore[indicates]=beta_combine
  ###########################################
  
  
  
  result_compare[[iter]] = cbind(FULL=coef(full_md),D_OPT=coef(dopt_train),LEV=coef(lev_train),UNIF=coef(rand_train),SPC=beta_restore)
  (nonzeros=unique(c(which(result_compare[[iter]][,1]!=0),which(result_compare[[iter]][,2]!=0),which(result_compare[[iter]][,3]!=0),which(result_compare[[iter]][,4]!=0))))
  out[[iter]] = as.data.frame(as.matrix(result_compare[[iter]][sort(nonzeros),]))
  names(out[[iter]]) = c("FULL","D_OPT","LEV","UNIF","SPC")
  
  #testing on new data
  (mse_fullpred = predmse(test1,full_md))
  (mse_dopt = predmse(test1,dopt_train))
  (mse_lev = predmse(test1,lev_train))
  (mse_unif = predmse(test1,rand_train))
   ( mse_spc = sum((as.matrix(cbind2(1,test1[,-281]))%*% as.matrix(beta_restore) - test1[,281])^2)/nrow(test1))

  
}
out
saveRDS(out,"bootstrap_real1000.rds")
saveRDS(result_compare,"bootstrap_real1000_full.rds")

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


temp = readRDS("result_compare_1000.rds")

