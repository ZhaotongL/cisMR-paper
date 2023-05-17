.libPaths(c(.libPaths(),"/panfs/jay/groups/20/panwei/lin00374/saonli/R/x86_64-pc-linux-gnu-library/3.6"))
library(MR.LDP)
library(data.table)
library(MendelianRandomization)
library(dplyr)
library(snpStats)
library(gsmr)
source('/home/panwei/lin00374/cisMR/LDA_Egger_simu/redo0414/cisMR_cML.R')
source('/home/panwei/lin00374/cisMR/LDAMR.R')


array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
job_id = as.character(Sys.getenv('SLURM_JOB_NAME'))

args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
      if(array_id==1){write.table(args[[i]],paste0('params',job_id,'.txt'),append=T,quote=F,row.names=F,col.names=F)}
    }
}



ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

LDA_Egger_simulation <- function(seed,h2e,h2y,gamma,N1,N2,rho,J,K1,K2,tau,p_causal){
  set.seed(seed)
  LD_mat = ar1_cor(J,rho)
  L = chol(LD_mat)
  beta_E = rep(0,J)
  causal.index = sample(1:J,J*p_causal)
  beta_E[causal.index] = rnorm(n=length(causal.index),mean=0,sd=1)
  sigma_E = sqrt(h2e/(t(beta_E) %*% LD_mat %*% beta_E))
  beta_E = c(sigma_E) * beta_E
  theta = rep(0,J)

  if(K1+K2>0){
    theta[c(sample(causal.index,K1),sample(setdiff(1:J,causal.index),K2))] = tau + rnorm(n=K1+K2,mean=0,sd=1)
    sigma = sqrt(h2y/(t(theta) %*% LD_mat %*% theta))
    theta = c(sigma) * theta
  }
  beta_G = theta + gamma*beta_E
  sx = rep(sqrt((1-h2e)/N1),J)
  sy = rep(sqrt((1-h2e*gamma^2-h2y)/N2),J)
  Sig_exp = (LD_mat*crossprod(t(sx))) 
  Sig_out = (LD_mat*crossprod(t(sy)))
  beta_E_hat = LD_mat %*% beta_E +  t(L) %*% rnorm(J,mean=0,sd=sx)
  beta_G_hat = LD_mat %*% beta_G +  t(L) %*% rnorm(J,mean=0,sd=sy)
  return(list(LD_mat=LD_mat,betax=beta_E_hat,sdx=sx,betay=beta_G_hat,sdy=sy,
              Sig_exp=Sig_exp,Sig_out=Sig_out,theta=theta,beta_E=beta_E,beta_G=beta_G, theta = theta,
              exp_IV = sort(causal.index)
  ))  
}

#N1 = 10000
#N2 = 80000
#J = 20
#p = 1
#K = 0
#byx = 0.2
#rho = 0.4

for(j in (20*(array_id-1)):(20*array_id-1)+1){
    s = LDA_Egger_simulation(seed=j,J=J,p_causal=p,K1=K1,K2=K2,N1 = N1,N2=N2,gamma = byx,rho=rho,h2e=h2x,h2y=h2y,tau=tau)
    b_exp = s$betax; b_out = s$betay; se_exp = s$sdx; se_out = s$sdy
    exp_IV = s$exp_IV
    R = s$LD_mat
    LD_inv = solve(s$LD_mat)
    b_exp_cond = LD_inv %*% b_exp; b_out_cond = LD_inv %*% b_out
    Sig_exp = LD_inv %*% (s$LD_mat * crossprod(t(se_exp))) %*% LD_inv;
    Sig_out = LD_inv %*% (s$LD_mat * crossprod(t(se_out))) %*% LD_inv;
    Sig_exp_inv = solve(Sig_exp); Sig_out_inv = solve(Sig_out); 
    stick0 = b_out_cond


   ciscML_res = cismr_cML_DP(b_exp_cond,b_out_cond,Sig_exp_inv,Sig_out_inv,maxit=200,n = N1,random_start = 5,random_start_pert = 5,method='naive')
   ciscML0_b = ciscML_res$BIC_theta;ciscML0_se=ciscML_res$BIC_se; ciscML0_pval=ciscML_res$BIC_p
   ciscML0_TP = length(intersect(ciscML_res$BIC_invalid,which(s$theta!=0)))
   ciscML0_FP = length(setdiff(ciscML_res$BIC_invalid,which(s$theta!=0)))
   ciscMLDP0_b = ciscML_res$BIC_DP_theta; ciscMLDP0_se=ciscML_res$BIC_DP_se; ciscMLDP0_pval=ciscML_res$BIC_DP_p
#
   LDAegger_res = LDA.MREgger(X=b_exp_cond,Y=b_out_cond,W=Sig_out_inv)
   LDAegger_b = LDAegger_res[2,1]; LDAegger_pval = LDAegger_res[2,4]; LDAegger_se = sqrt(LDAegger_res[2,6])

    mrinput = mr_input(bx=as.vector(b_exp), bxse=as.vector(se_exp), by=as.vector(b_out), byse=as.vector(se_out), correlation=R)
   ivw_res = mr_ivw(mrinput,correl=TRUE)
   ivw_b = ivw_res@Estimate; ivw_se = ivw_res@StdError; ivw_pval=ivw_res@Pvalue
    egger_res = mr_egger(mrinput,correl=TRUE)
    egger_b = egger_res@Estimate; egger_se = egger_res@StdError.Est; egger_pval = egger_res@Pvalue.Est



    exp_IV = s$exp_IV
    b_exp = s$betax[exp_IV]; b_out = s$betay[exp_IV]; se_exp = s$sdx[exp_IV]; se_out = s$sdy[exp_IV]
    R = s$LD_mat[exp_IV,exp_IV]
    LD_inv = solve(R)
    b_exp_cond = LD_inv %*% b_exp; b_out_cond = LD_inv %*% b_out
    Sig_exp = LD_inv %*% (R * crossprod(t(se_exp))) %*% LD_inv;
    Sig_out = LD_inv %*% (R * crossprod(t(se_out))) %*% LD_inv;
    Sig_exp_inv = solve(Sig_exp); Sig_out_inv = solve(Sig_out);
    stick0 = b_out_cond


   LDAeggerX_res = LDA.MREgger(X=b_exp_cond,Y=b_out_cond,W=Sig_out_inv)
   LDAeggerX_b = LDAeggerX_res[2,1]; LDAeggerX_pval = LDAeggerX_res[2,4]; LDAeggerX_se = sqrt(LDAeggerX_res[2,6])
#
#    LDAivwX_res = LDA.MRIVW(X=b_exp_cond,Y=b_out_cond,W=Sig_out_inv)
#    LDAivwX_b = LDAivwX_res[1,1]; LDAivwX_pval = LDAivwX_res[1,5]; LDAivwX_se = sqrt(LDAivwX_res[1,2])
#    LDAivwXt_pval = LDAivwX_res[1,4]
#
    mrinput = mr_input(bx=as.vector(b_exp), bxse=as.vector(se_exp), by=as.vector(b_out), byse=as.vector(se_out), correlation=R)
   ivwX_res = mr_ivw(mrinput,correl=TRUE)
   ivwX_b = ivwX_res@Estimate; ivwX_se = ivwX_res@StdError; ivwX_pval=ivwX_res@Pvalue
    eggerX_res = mr_egger(mrinput,correl=TRUE)
    eggerX_b = eggerX_res@Estimate; eggerX_se = eggerX_res@StdError.Est; eggerX_pval = eggerX_res@Pvalue.Est
#
   ciscMLX_res = cismr_cML_DP(b_exp_cond,b_out_cond,Sig_exp_inv,Sig_out_inv,maxit=200,n = N1,random_start = 5,random_start_pert = 5,method='naive')
   ciscMLX0_b = ciscMLX_res$BIC_theta; ciscMLX0_se=ciscMLX_res$BIC_se; ciscMLX0_pval=ciscMLX_res$BIC_p
   ciscMLXDP0_b = ciscMLX_res$BIC_DP_theta; ciscMLXDP0_se=ciscMLX_res$BIC_DP_se; ciscMLXDP0_pval=ciscMLX_res$BIC_DP_p
#


    
    var_name =  c(ls(pattern='_b$'),ls(pattern='_se$'),ls(pattern='_pval$'),ls(pattern='_TP$'),ls(pattern='_FP$'))
    for(vn in var_name){
        t1 = paste0("write.table(t(",vn,"),'",vn,"_",job_id,"_",array_id,".txt',quote=F,row.names=F,append=T,col.names=F)")
        eval(parse(text=t1))
    }
}
