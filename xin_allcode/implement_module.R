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
#Function Description: this function returns list of regression parameter/time consumption and so on
#Dependencies: 
#Input: 
#Output:
#######################################################################
reg_module<- function(#data will be passed from global environment
                      method,#approach including "all" "cor" "mod" "full"
                      r=cor_prop,
                      dist_resp = "gaussian", #applicable for logistic regression
                      seed_l=seed6,#seed for leverage sampling
                      ...#arguments needed by specific method
){
  
  if(! method %in% c("cor","all","mod","full","rand","lev")) {
    print("error!")
    return(NA)} 
  if(method == "cor"){row_select = cor_var_max}
  if(method == "all"){row_select = all_vars_max}
  if(method == "mod"){row_select = mod_var_max}
  if(method == "lev"){row_select = lev_sample}
  if(method == "rand"){
    select_time<- system.time(rand_idx <- sample(1:dim(x)[1],subsize)) #this x need to have index in the first column
    rand_x <- as.matrix(x_mat[rand_idx,])
    rand_y <-  as.matrix(y[rand_idx])
    lasso_time<- system.time(md_cv <- cv.glmnet(rand_x,rand_y,alpha=islasso,family = dist_resp))
    coef_md_cv = coef(md_cv)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    n_selected = subsize
  }
  #IBOSS subset and its lasso regression
  if(method %in% c("cor","all","mod","lev")){
    select_time<- system.time(informative_idx <- row_select(x,y,subsize,prop = r,seed=seed_l))
    sub_x <- as.matrix(x_mat[informative_idx,])
    sub_y = as.matrix(y[informative_idx])
    (n_selected = length(informative_idx))#to make sure random sample and informative subset have same size
    lasso_time<- system.time(md_cv <- cv.glmnet(sub_x,sub_y,alpha=islasso,family = dist_resp))
    coef_md_cv = coef(md_cv)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    (n_selected = length(informative_idx))
  }
  if(method == "full"){
    select_time = c(0,0,0)
    lasso_time<- system.time(md_cv <- cv.glmnet(x_mat,y,alpha=islasso,family = dist_resp))
    coef_md_cv = coef(md_cv)
    mse.min <-md_cv$cvm[md_cv$lambda == md_cv$lambda.min]
    n_selected = n
  }
    return(list(method=method,proprotion=r,coef = coef_md_cv,mse=mse.min,slct_time = select_time[3],lasso_time =lasso_time[3],n_slct=n_selected))
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
  MPE = mse(y_hat,y_test)
  report_row = data.frame(s = num_catch(coef_l),
             TPR = sensitivity(coef_l),
             SPC = specificity(coef_l),
             PPV = precision(coef_l),
             ACC = accuracy(coef_l),
             MPE = MPE,
             MSE = reg_list$mse,
             n_slct= reg_list$n_slct,
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