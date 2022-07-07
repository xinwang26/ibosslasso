

split_n_conquer <- function(num_slice=5,vote_rate=0.4,dist_resp = "gaussian", para_set = F,...){

  omega = floor(vote_rate*num_slice)
  #slice data
  #######################################################################
  nk = floor(n/num_slice)  #num_slice is just k in the paper, name k is probably occupied, incase override use num_slice
  n_selected = nk
  nk_end = n - (num_slice-1)*nk #make sure the last share cover the remain
  registerDoParallel(num_slice)
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
              ,proprotion=r
              ,coef = beta_restore
              ,mse = mse.min  #avg of the 5
              ,slct_time = comb_time[3]  #combination time for spc
              ,lasso_time = lasso_time[3]
              ,n_slct = n_selected #nk
              ,cvm_rch_min = cvm_rch_min 
              ) 
         )
}
