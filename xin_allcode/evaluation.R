#######################################################################
#Author: Xin Wang
#Initial Date: Fri Jun 10 15:17:51 2016
#Program Description: this program contains all function computing the four columns to evaluate performanc
#Completed: 
#Dependencies: 
#Main Reference: 
#Input: 
#Output:
#######################################################################

num_catch <-function(coef_l){return(sum(coef_l[-1]!=0))} #this function calculate number of non-zero estimates
sensitivity <-function(coef_l,true_b=true_beta,p_t=p_true){sum(coef_l*true_b!=0)/p_t}
specificity <-function(coef_l,true_b=true_beta,p_t=p_true,pb=p){sum(coef_l==0&true_b==0)/(pb+1-p_t)}
precision <- function(coef_l,true_b=true_beta,p_t=p_true){sum(coef_l*true_b!=0)/sum(coef_l!=0)}
accuracy <- function(coef_l,true_b=true_beta,p_t=p_true,pb=p+1){( sum(coef_l*true_b!=0) + sum(coef_l==0&true_b==0 ))/pb }
sign_consistency <-function(coef_l,true_b=true_beta,p_t=p_true,pb=p+1){(sum(coef_l*true_b!=0)+sum(coef_l==0&true_b==0))/p}
train_mse<-function(x_x,x_y,coef_l){
  y_hat = x_x%*%matrix(coef_l[-1],ncol=1)+ matrix(rep(1,dim(x_x)[1]),ncol =1)*coef_l[1]
  train_mse = mse(y_hat,x_y)
  return(train_mse)
}


test_mse<-function(coef_l,size = 100,seed5 = 12332454,data_generator = data_gen_n,p_full =p ,beta_t = beta,sig = sigma1){
  x_test <- data_generator(seed5,n_full = size)
  x_test <- as.matrix(x_test[,paste0("V",1:p_full),with=F])
  y_test <- x_test%*%beta_t+rnorm(size,0,sig) 
  (a = train_mse(x_test,y_test,coef_l))
  return(a)
}

sse_beta <-function(coef_l,beta = true_beta){return(sum((coef_l-beta)^2))}

repo_combine<-function(coef_l){
  return(
  data.frame(s = num_catch(coef_l),
        TPR = sensitivity(coef_l),
        SPC = specificity(coef_l),
        PPV = precision(coef_l), 
        ACC = accuracy(coef_l),
        MSE = test_mse(coef_l,test_size,seed5,data_gen)),
        sign_consi = sign_consistency(coef_l)
  )
}
#positive_catch <-function(coef_l,true_beta = true_beta){return(sum(true_beta*coef_l!=0))} #this function calculate the number of that catching the true non-zero coefficient
#sse<-function(coef_l,true_b = true_beta){return(sum((coef_l-true_b)^2))}
# sensitive_catch <- function(coef_mat,true_b = true_beta,p_t=p_true){#this function find out first column that catching equal number variable to the true, calculate its sensitivity
#   (num_of_catch = apply(coef_mat,2,num_catch))
#   idx = which(num_of_catch>=p_t)[1]
#   ests = coef_mat[,idx]
#   num_catch = sum(ests!=0)
#   num_true_catch = sum(ests*true_b!=0)
#   (rate =num_true_catch/p_true)
#   return(list(idx = idx,num_catch = num_catch, num_true_catch = num_true_catch,sensitivity = rate,est = ests))
# }
# 
# specific_catch <- function(coef_mat,true_beta1 = true_beta, p_t=p_true){#this function find the first column that catching all the true variables and calculate its specificity
#   (num_to_catch = apply(coef_mat,2,positive_catch,true_beta=true_beta1))
#   idx = which(num_to_catch>=p_t)[1]
#   ests = coef_mat[,idx]
#   num_catch = sum(ests!=0)
#   num_true_catch = sum(ests*true_beta!=0)
#   (rate =(p-num_catch) /(p-sum(true_beta1!=0)))
#   return(list(idx = idx,est = ests,num_catch = num_catch, num_true_catch = num_true_catch,specificity = rate))
# }

