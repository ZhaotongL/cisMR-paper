.libPaths(c(.libPaths(),"/panfs/jay/groups/20/panwei/lin00374/saonli/R/x86_64-pc-linux-gnu-library/3.6"))
source('/home/panwei/lin00374/cisMR/LDA_Egger_simu/redo0414/cisMR_cML.R')
array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


Res_df = NULL
for(f in list.files(pattern='.RData')[(30*(array_id-1)):(30*array_id-1)+1]){
    if(is.na(f)){next;}
    s = unlist(strsplit(f,'.RData'))
    load(f)
    b_exp_cond=mr_dat$b_exp
    b_out_cond=mr_dat$b_out
    mr_dat$out_df$se_cor[which(is.na(mr_dat$out_df$se_cor))] = sqrt(1/mr_dat$out_df$N[1])
    Sig_exp1 = solve(mr_dat$LD_mat) %*% (mr_dat$exp_df$se_cor %o% mr_dat$exp_df$se_cor * mr_dat$LD_mat) %*% solve(mr_dat$LD_mat)
    Sig_out1 = solve(mr_dat$LD_mat) %*% (mr_dat$out_df$se_cor %o% mr_dat$out_df$se_cor * mr_dat$LD_mat) %*% solve(mr_dat$LD_mat)
    #Sig_exp1 = mr_dat$Sig_exp
    #Sig_out1 = mr_dat$Sig_out
    Sig_exp_inv=solve(Sig_exp1)
    Sig_out_inv=solve(Sig_out1)

    stick0 = mr_dat$b_out # for the sake of my code
    set.seed(12345)
    ciscML_res = try(cismr_cML_DP(b_exp=b_exp_cond,b_out=b_out_cond,Sig_exp_inv=Sig_exp_inv,Sig_out_inv=Sig_out_inv,maxit=200,n = mr_dat$N1,random_start = 5,min_theta_range=-0.1,max_theta_range=0.1,num_pert=100,random_start_pert=5))
    if(class(ciscML_res)=='try-error'){print(paste0('Error: ',s));next;}
    write.table(t(c(s,length(stick0),ciscML_res$BIC_DP_theta,ciscML_res$BIC_DP_se,ciscML_res$BIC_DP_p,length(ciscML_res$BIC_invalid))),paste0('cisMR_res.txt',array_id),quote=F,row.names=F,col.names=F,append=T)
    write.table(t(c(s,ciscML_res$BIC_theta,ciscML_res$BIC_se,ciscML_res$BIC_p)),paste0('cisMRBIC_res.txt',array_id),quote=F,row.names=F,col.names=F,append=T)
}

