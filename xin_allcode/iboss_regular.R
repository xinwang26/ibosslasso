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
setwd("/Users/xinwang/Google Drive/research reading/2016Spring/BIGDATA/code")
library(data.table)
library(Metrics)

library(lars)
library(glmnet)
source("data_Gen.R",echo = T,max.deparse.length = 10000)
#this contains functions for the ranking and selecting
source("data_table_funcs_small.R",echo = T,max.deparse.length = 10000)
#library(irlba)


datetime = substr(Sys.time(),6,16)
n = 5e3
p = 5e2
subsize = 1e3
test_size = 100
cor_prop = 0.20
n_true = sqrt(p)+1
rounds = 100
set.seed(3)
seed1_list = sample(1:99999,rounds)
seed2_list = sample(99999:199999,rounds)
seed3_list = sample(200000:299999,rounds)
seed4_list = sample(300000:399999,rounds)
seed4_list = sample(400000:499999,rounds)
seed5_list = sample(500000:599999,rounds)
(scale_beta  = 1.5*sqrt(log(p)/subsize))#necessary for data generation

X_distribution = "normal"
data_gen = data_gen_n
row_selectA = all_vars_max #equally consider all variables / cor_var_max/ mod_var_max/cor_var_min/cor_to_max/all_vars_max
row_selectC = cor_var_max
sigma1 = 1 #scale of error


source("evaluation.R") # the evaluation function test-mse will rely on the distribution of predictors

#sens = data.frame()
#spec = data.frame()
#prec = data.frame()
#accr = data.frame()
#mses = data.frame()
#nums = data.frame()
#sseb = data.frame()
#temses = data.frame()
#trmses = data.frame()
ptm = proc.time()
full_repo = data.frame()
subA_repo = data.frame()
subC_repo = data.frame()
rand_repo = data.frame()

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
 
  #generate the full data set
  system.time(x_mat <- as.matrix(x[,paste0("V",1:p),with=F]))
  set.seed(seed3)
  y <-  as.matrix(x_mat)%*%beta+rnorm(n,0,sigma1)
  #perform lasso regression on full data
  md_time<- system.time(md_cv <- cv.glmnet(x_mat,y,alpha=1,family = "gaussian"))
  coef_md_cv = coef(md_cv)
  
  #IBOSS subset and its lasso regression
  subA_select_time<- system.time(informativeA_idx <- row_selectA(x,y,subsize))
  subA_x <- as.matrix(x_mat[informativeA_idx,])
  subA_y = as.matrix(y[informativeA_idx])
  (n_selectedA = length(informativeA_idx))#to make sure random sample and informative subset have same size
  subA_time<- system.time(subA_md_cv <- cv.glmnet(subA_x,subA_y,alpha=1,family = "gaussian"))
  coef_subAmd_cv = coef(subA_md_cv)
  
  
  #IBOSS subset and its lasso regression
  subC_select_time<- system.time(informativeC_idx <- row_selectC(x,y,subsize,prop = cor_prop))
  subC_x <- as.matrix(x_mat[informativeC_idx,])
  subC_y = as.matrix(y[informativeC_idx])
  (n_selectedC = length(informativeC_idx))#to make sure random sample and informative subset have same size
  subC_time<- system.time(subC_md_cv <- cv.glmnet(subC_x,subC_y,alpha=1,family = "gaussian"))
  coef_subCmd_cv = coef(subC_md_cv)
  
  
  #random subset and its lasso regression
  set.seed(seed4)
  rand_select_time<- system.time(rand_idx <- sample(1:dim(x)[1],n_selectedC)) #this x need to have index in the first column
  rand_x <- as.matrix(x_mat[rand_idx,])
  rand_y <-  as.matrix(y[rand_idx])
  rand_time<- system.time(rand_md_cv <- cv.glmnet(rand_x,rand_y,alpha=1,family = "gaussian"))
  coef_randmd_cv = coef(rand_md_cv)
  
  
  full_repo <- rbind(full_repo, 
                     c(repo_combine(coef_md_cv),lasso_time = md_time[3],slct_time.elapsed = 0))  #six measurements and time consumned
  subA_repo <- rbind(subA_repo,
                     c(repo_combine(coef_subAmd_cv),lasso_time = subA_time[3],slctA_time = subA_select_time[3]))
  subC_repo <- rbind(subC_repo,
                    c(repo_combine(coef_subCmd_cv),lasso_time = subC_time[3],slctC_time = subC_select_time[3]))
  rand_repo <- rbind(rand_repo,
                     c(repo_combine(coef_randmd_cv),lasso_time = rand_time[3],slct_time = rand_select_time[3]))
  #implemented by lars:
  #system.time(cv_sub <- cv.lars(sub_x,sub_y,type= "lasso"))
  #system.time(sub_md_cv <- lars(sub_x,sub_y,type= "lasso"))
  #(l_sub = cv_sub$index[which.min(cv_sub$cv)])
  #coef_submd_cv = predict(sub_md_cv,sub_x,s=l_sub,type="coefficient",mode="fraction")$coef
  
  #system.time(cv_rand <- cv.lars(rand_x,rand_y,type= "lasso"))
  #system.time(rand_md_cv <- lars(rand_x,rand_y,type= "lasso"))
  #(l_rand = cv_rand$index[which.min(cv_rand$cv)])
  #coef_randmd_cv = predict(rand_md_cv,rand_x,s=l_rand,type="coefficient",mode="fraction")$coef
  
  #system.time(cv_md <- cv.lars(x_mat,y,type= "lasso"))
  #system.time(md_cv <- lars(x_mat,y,type= "lasso"))
  #(l_md = cv_md$index[which.min(cv_md$cv)])
  #coef_md_cv = predict(md_cv,x_mat,s=l_md,type="coefficient",mode="fraction")$coef
  

    
  #(nums = rbind(nums,c(num_catch(coef_md_cv),
  #                     num_catch(coef_submd_cv ),
   #                    num_catch(coef_randmd_cv))))
  #(sens = rbind(sens,c(sensitivity(coef_md_cv),
   #                    sensitivity(coef_submd_cv ),
    #                   sensitivity(coef_randmd_cv))))
  #(spec = rbind(spec,c(specificity(coef_md_cv),
   #                    specificity(coef_submd_cv ),
    #                   specificity(coef_randmd_cv))))
  #(accr = rbind(accr,c(accuracy(coef_md_cv),
   #                    accuracy(coef_submd_cv ),
    #                   accuracy(coef_randmd_cv))))
  #(prec = rbind(prec,c(precision(coef_md_cv),
   #                    precision(coef_submd_cv ),
    #                   precision(coef_randmd_cv))))
  #(temses = rbind(temses,c(test_mse(coef_md_cv,test_size,seed5,data_generator = data_gen),
   #                        test_mse(coef_submd_cv,test_size,seed5,data_generator = data_gen),
    #                       test_mse(coef_randmd_cv,test_size,seed5,data_generator = data_gen))))
  #(trmses = rbind(trmses,c(train_mse(x_mat,y,md_cv),
   #                        train_mse(x_mat,y,sub_md_cv),
    #                       train_mse(x_mat,y,rand_md_cv))))
  #(sseb = rbind(sseb,c(sse_beta(coef_md_cv),
   #                    sse_beta(coef_submd_cv ),
    #                   sse_beta(coef_randmd_cv))))
  
}





#write.csv(sens,file = paste0("sensi",datetime,".csv"))
#write.csv(spec,file = paste0("spec",datetime,".csv"))
#write.csv(temses,file = paste0("trmses",datetime,".csv"))
#write.csv(trmses,file = paste0("trmses",datetime,".csv"))
#write.csv(sseb,file = paste0("mses",datetime,".csv"))
sink(file=paste0("reg_R_result_table",datetime,".txt"))
cat("All the following result are implemented by R package glmnet")
full_mean = round(apply(full_repo,2,mean),4)
full_sd = round(apply(full_repo, 2, sd),4)
subA_mean = round(apply(subA_repo, 2, mean),4)
subA_sd = round(apply(subA_repo, 2, sd),4)
subC_mean = round(apply(subC_repo, 2, mean),4)
subC_sd = round(apply(subC_repo, 2, sd),4)
rand_mean = round(apply(rand_repo, 2, mean),4)
rand_sd = round(apply(rand_repo, 2, sd),4)

output = c()
loop_all = c("full","subC","subA","rand") 
for(case in loop_all){
  assign(paste0(case,"_tab"),c())
  for(j in 1:length(get(paste0(case,"_mean"))) ){
    temp = eval(parse(text = paste0("c(",case,"_mean[",j,"],","sep='&(',",case,"_sd[",j,"]",",sep=')&'",")")))
    assign(paste0(case,"_tab"), 
           eval(parse(text=paste0("c(",paste0(case,"_tab,temp)") )))
           )
  }
  output = rbind(output,get(paste0(case,"_tab")))
}
rownames(output)<-loop_all

conditions = data.frame(X_distribution,rounds,n,p,n_true,scale_beta,subsize,n_selectedA,cor_prop,n_selectedC)
print(conditions)
print(output)
sink()
#output to csv
write.table("",sep=",",file= paste0("reg_R_result_table",datetime,".csv"))
write.table(conditions,sep=",",file = paste0("reg_R_result_table",datetime,".csv"),append=T,col.names = T)
write.table(output,sep = ",",file = paste0("reg_R_result_table",datetime,".csv"),append = T)
