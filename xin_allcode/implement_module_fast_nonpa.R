#######################################################################
#Author: Xin Wang
#Initial Date: Thu Sep 22 20:17:27 2016
#Program Description: this script contains functions needed for the one full implement 
#Completed: 
#Dependencies: data_table_funcs_small.R, require will be put in the main loop script
#Main Reference: 
#Input: 
#Output:
#######################################################################
#######################################################################
#Function Description: 
#Dependencies: 
#Input: 
#Output:regression with only correlated columns
#######################################################################
sis <- function(#data will be passed from global environment
  method,#approach including "all" "cor" "mod" "full"
  r=0.05, #ratio of selected for SIS, vote ratio if split and conquer
  dist_resp = "gaussian", #applicable for logistic regression
  seed_l=seed6,#seed for leverage sampling
  ...#arguments needed by specific method
)
{
  if(! method %in% c("cor",'clev','calev','full','rand')) {  #subdata approaches
    print("error!")
    return(NA)} 
  if(method == "cor"){row_select = cor_var_max_sis}
  if(method == "clev"){row_select = cor_lev_sample_sis}
  if(method == "calev"){row_select = cor_alev_sample_sis}
  if(method == "rand"){
    (to_select_size = floor(dim(x)[2]*r))
    corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
    corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
    select_time = 0
    slct_vars = to_select
    set.seed(seed_l)
    select_time<- system.time(rand_idx <- sample(1:dim(x)[1],subsize)) #this x need to have index in the first column
    select_time = select_time[3]
    rand_x <- as.matrix(x_mat[rand_idx,slct_vars])
    rand_y <-  as.matrix(y[rand_idx])
    lasso_time<- system.time(md_cv <- cv.glmnet(rand_x,rand_y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    temp_coef = coef(md_cv)
    buffer =  as(matrix(rep(0,p),ncol=1), "dgCMatrix")
    buffer[slct_vars,1] =  temp_coef[-1]
    coef_md_cv = rbind(temp_coef[1],buffer)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    (n_selected = length(rand_idx))
  }
  if(method == "full"){
    (to_select_size = floor(dim(x)[2]*r))
    corcalc_time <- system.time(cor_resp <- unlist(lapply(x[,-dim(x)[2],with=F],cor,y=y)))
    corsort_time<- system.time(to_select <- sort(abs(cor_resp),decreasing = T,index.return = T)$ix[1:to_select_size])
    select_time = 0
    slct_vars = to_select
    sub_x <- as.matrix(x_mat[,slct_vars])
    set.seed(seed_l)
    lasso_time<- system.time(md_cv <- cv.glmnet(sub_x,y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    temp_coef = coef(md_cv)
    buffer =  as(matrix(rep(0,p),ncol=1), "dgCMatrix")
    buffer[slct_vars,1] =  temp_coef[-1]
    coef_md_cv = rbind(temp_coef[1],buffer)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    (n_selected = n)
    }
  if(method %in% c("cor","all","mod","lev","inno",'alev','clev','calev')){
    rows_n_cols <- row_select(x,y,subsize,prop = r,seed=seed_l)
    select_time <- rows_n_cols$select_time
    informative_idx = rows_n_cols$rows
    slct_vars = rows_n_cols$cols
    sub_x <- as.matrix(x_mat[informative_idx,slct_vars])
    sub_y = as.matrix(y[informative_idx])
    set.seed(seed_l)
    lasso_time<- system.time(md_cv <- cv.glmnet(sub_x,sub_y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    temp_coef = coef(md_cv)
    buffer =  as(matrix(rep(0,p),ncol=1), "dgCMatrix")
    buffer[slct_vars,1] =  temp_coef[-1]
    coef_md_cv = rbind(temp_coef[1],buffer)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    (n_selected = length(informative_idx))
  }
  return(list(method=method,proprotion=r,coef = coef_md_cv,mse=mse.min,slct_time = select_time,lasso_time =lasso_time[3],n_slct=n_selected,cvm_rch_min = (md_cv$lambda.min!=md_cv$lambda[length(md_cv$lambda)])))
}


#######################################################################
#Function Description: split and conquer approach regression and result
#Dependencies: 
#Input: global variables, slice number and voting rate, parallel setting
#Output:
#######################################################################
split_n_conquer <- function(num_slice=5,vote_rate=0.4,dist_resp = "gaussian", para_set = F,...){
  method = "spc"
  omega = floor(vote_rate*num_slice)
  #slice data
  #######################################################################
  nk = floor(n/num_slice)  #num_slice is just k in the paper, name k is probably occupied, incase override use num_slice
  n_selected = nk
  nk_end = n - (num_slice-1)*nk #make sure the last share cover the remain
  registerDoParallel(1)
  #parallel compute
  #######################################################################
  #since the data itself is totally randomlly generated, no need to sort and do random subset, simple slicing is enough
  ptm_start = proc.time()
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
    md_cv <- cv.glmnet(slice_x,slice_y,parallel = para_set,alpha=islasso,lambda.min.ratio = lam_min_rate,family = dist_resp)
    cvm_rch_min = (md_cv$lambda.min!=md_cv$lambda[length(md_cv$lambda)])
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    betak=list(betahat=coef(md_cv),mse.min=mse.min,cvm_rch_min = 1*cvm_rch_min)
  }
  lasso_time = proc.time()-ptm_start
  coef_list = lapply(split_md, `[[`, 1)
  split_beta = do.call(cbind, unlist(coef_list, recursive=TRUE))
  mse.min =  mean(unlist(lapply(split_md, `[[`, 2), recursive=TRUE))
  cvm_rch_min = min(unlist(lapply(split_md, `[[`, 3), recursive=TRUE))
  #combine
  #######################################################################
  theta = list()
  nonzero_ind= as.matrix(1*(split_beta!=0))
  tic = proc.time()
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
    Xk = cbind(rep(1,n_k), x_mat[slice_start:slice_end,] )[,indicates]
    system.time(theta[[i]] <- Xk %*% split_beta[indicates,i] )
    system.time(component[[i]] <- E%*%t(Xk)%*%Xk%*%E)
    system.time(betai[[i]] <- component[[i]]%*%split_beta[indicates,i])
  }
  system.time(beta_combine <- E %*% ginv(as.matrix(Reduce(`+`, component))) %*% Reduce(`+`, betai))
  comb_time = proc.time()-tic
  beta_restore = as(matrix(rep(0,p+1),ncol=1),"dgCMatrix")
  beta_restore[indicates]=beta_combine
  return(list(method=method
              ,proprotion=vote_rate
              ,coef = beta_restore
              ,mse = mse.min  #avg of the 5
              ,slct_time = comb_time[3]  #combination time for spc
              ,lasso_time = lasso_time[3]
              ,n_slct = n_selected #nk
              ,cvm_rch_min = cvm_rch_min 
  ) 
  )
}
#######################################################################
#Function Description: this function returns list of regression parameter/time consumption and so on
#Dependencies: 
#Input: 
#Output:
#######################################################################
reg_module<- function(#data will be passed from global environment
                      method,#approach including "all" "cor" "mod" "full"
                      r=cor_prop, #ratio of selected for SIS, vote ratio if split and conquer
                      dist_resp = "gaussian", #applicable for logistic regression
                      seed_l=seed6,#seed for leverage sampling
                      ...#arguments needed by specific method
){
  
  if(! method %in% c("cor","all","mod","full","rand","lev","alev",'clev','calev',"inno","spc")) {  #subdata approaches
    print("error!")
    return(NA)} 
  if(method == "cor"){row_select = cor_var_max}
  if(method == "all"){row_select = all_vars_max}
  if(method == "mod"){row_select = mod_var_max}
  if(method == "lev"){row_select = lev_sample}
  if(method == "alev"){row_select = approxlev_sample}
  if(method == "clev"){row_select = cor_lev_sample}
  if(method == "calev"){row_select = cor_alev_sample}
  if(method == "rand"){
    set.seed(seed_l)
    select_time<- system.time(rand_idx <- sample(1:dim(x)[1],subsize)) #this x need to have index in the first column
    rand_x <- as.matrix(x_mat[rand_idx,])
    rand_y <-  as.matrix(y[rand_idx])
    lasso_time<- system.time(md_cv <- cv.glmnet(rand_x,rand_y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    coef_md_cv = coef(md_cv)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    n_selected = subsize
  }#if(method == "inno"){row_select = two_step_range;dist_resp = "binomial"}
  #IBOSS subset and its lasso regression
  if(method %in% c("cor","all","mod","lev","inno",'alev','clev','calev')){
    select_time<- system.time(informative_idx <- row_select(x,y,subsize,prop = r,seed=seed_l))
    sub_x <- as.matrix(x_mat[informative_idx,])
    sub_y = as.matrix(y[informative_idx])
    set.seed(seed_l)
    lasso_time<- system.time(md_cv <- cv.glmnet(sub_x,sub_y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    coef_md_cv = coef(md_cv)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    (n_selected = length(informative_idx))
  }
  if(method == "full"){
    select_time = c(0,0,0)
    lasso_time<- system.time(md_cv <- cv.glmnet(x_mat,y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    coef_md_cv = coef(md_cv)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    n_selected = n
  }
  if(method == "spc"){
    return(split_n_conquer(num_slice=5,vote_rate=r,dist_resp = "gaussian", para_set = F))
    }
    return(list(method=method,proprotion=r,coef = coef_md_cv,mse=mse.min,slct_time = select_time[3],lasso_time =lasso_time[3],n_slct=n_selected,cvm_rch_min = (md_cv$lambda.min!=md_cv$lambda[length(md_cv$lambda)])))
}

#######################################################################
#Function Description: this function calculate the smallest number to keep in step of SIS so that the final SIS-IBOSS-LASSO model could include the true model
#Dependencies: p^* p, but x,y, need to be defined already in the global memory
#Input: 
#Output:
#######################################################################
s_find<-function(pstar,p,seed_l=seed4,dist_resp= "gaussian", ...)
{
  sensi = 0
  s = pstar
  while(sensi<1)
  {
    row_select = cor_var_max 
    select_time<- system.time(informative_idx <- row_select(x,y,subsize,prop = s/p,seed=seed_l))
    sub_x <- as.matrix(x_mat[informative_idx,])
    sub_y = as.matrix(y[informative_idx])
    set.seed(seed_l)
    lasso_time<- system.time(md_cv <- cv.glmnet(sub_x,sub_y,parallel = F,alpha=islasso,family = dist_resp,lambda.min.ratio = 0.0005))
    coef_md_cv = coef(md_cv)
    sensi <- sensitivity(coef_md_cv)
    s <- s+1
  }
  return(s)
}

#######################################################################
#Function Description: this function compute all the measurment of performance to be reported
#Dependencies: 
#Input: 
#Output:
#######################################################################

eval_module<- function(reg_list)
{
  (coef_l = reg_list$coef)
  y_hat = x_test%*%matrix(coef_l[-1],ncol=1)+ matrix(rep(1,dim(x_test)[1]),ncol =1)*coef_l[1]
  EPE = mse(y_hat,y_true)
  MPE = mse(y_hat,y_test)
  report_row = data.frame(s = num_catch(coef_l),
             TPR = sensitivity(coef_l),
             SPC = specificity(coef_l),
            # PPV = precision(coef_l),
            # ACC = accuracy(coef_l),
             MPE = MPE,
            EPE = EPE,
            # MSE = reg_list$mse,
             n_slct= reg_list$n_slct,
             sign_consi = sign_consistency(coef_l),
             slct_time = reg_list$slct_time,
             lasso_time = reg_list$lasso_time
             ,row.names = paste0(reg_list$method,reg_list$proprotion))
  return(report_row)
}

#######################################################################
#Function Description: this function get the looping result's mean and standard errors
#Dependencies: 
#Input: a data frame contains measuments from all loops of one specific setup
#Output: a df of one row containing all measurements' mean and sd
#######################################################################
report_module<-function(repo_df)
{
  repo_mean = round(apply(repo_df,2,mean),4)
  repo_sd = round(apply(repo_df,2,sd),4)
  report = c()
  name_char = names(repo_df)
  name_repo = c()
  for(j in 1:(length(repo_mean)-1) ){
    name_repo = cbind(name_repo,paste0(name_char[j],"_mean"),"__",paste0(name_char[j],"_sd"),"__")
    report = cbind(report,repo_mean[j],"&(",repo_sd[j],")&")
  }
  report = data.frame(cbind(report,repo_mean[j+1],"&(",repo_sd[j+1],")\\"))
  name_repo = cbind(name_repo,paste0(name_char[j+1],"_mean"),"__",paste0(name_char[j+1],"_sd"),"__")
  names(report) = name_repo
  return(report)
}

reportmean_module<-function(repo_df)
{
  repo_mean = round(apply(repo_df,2,mean),4)
  report = c()
  name_char = names(repo_df)
  name_repo = c()
  for(j in 1:(length(repo_mean)-1) ){
    name_repo = cbind(name_repo,paste0(name_char[j],"_mean"),"__")
    report = cbind(report,repo_mean[j],"&")
  }
  report = data.frame(cbind(report,repo_mean[j+1],"\\"))
  name_repo = cbind(name_repo,paste0(name_char[j+1],"_mean"),"__")
  names(report) = name_repo
  return(report)
}