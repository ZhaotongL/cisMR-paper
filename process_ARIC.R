.libPaths(c(.libPaths(),"/panfs/jay/groups/20/panwei/lin00374/saonli/R/x86_64-pc-linux-gnu-library/3.6"))
library(data.table)
library(dplyr)
source("/home/panwei/lin00374/TWAS_XSquare/Feb4_2021/allele_qc.R")


map = fread('~/EA/ARIC_pQTL.csv')
protein_list = names(which(table(map$SOMAmerID)>=3))
map = map %>% filter(SOMAmerID %in% protein_list)

array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


gwas = fread('~/cisMR/realdata/CAD/29212778-GCST005194-EFO_0000378-build37.f.tsv.gz')
gwas = gwas %>% select(SNP=variant_id, CHR=chromosome, A1=effect_allele, A2=other_allele, EAF=effect_allele_frequency, BETA=beta, SE=standard_error, P=p_value)
gwas$N = 547261
gwas$A1 = toupper(gwas$A1)
gwas$A2 = toupper(gwas$A2)

#a = list.files(path='./MRdat/',pattern='.RData')
for(i in (30*(array_id-1)):(30*array_id-1)+1){
    if(i>length(protein_list)){next;}
    s = protein_list[i]
    t = fread(paste0('~/EA/',s,'.PHENO1.glm.linear'))
    t$A0=t$REF
    t$A0[which(t$ALT!=t$A1 & t$A0==t$A1)] = t$ALT[which(t$ALT!=t$A1 & t$A0==t$A1)]
    chr = t$`#CHROM`[1]

    pqtl_df = data.frame(SNP=t$ID, A1=t$A1, A2=t$A0, freq=t$A1_FREQ, b=t$BETA, se=t$SE, p=t$P, N=t$OBS_CT)
    gwas_df = gwas %>% filter(SNP %in% pqtl_df$SNP)
    pqtl_df = pqtl_df %>% filter(SNP %in% gwas_df$SNP)
    gwas_pqtl_df = merge(gwas_df,pqtl_df,by='SNP') 
    

    remove_flip = allele.qc(gwas_pqtl_df$A1.x,gwas_pqtl_df$A2.x,
                          gwas_pqtl_df$A1.y,gwas_pqtl_df$A2.y)
    flip_snp = gwas_pqtl_df$SNP[remove_flip$flip]
    gwas_pqtl_df$BETA[which(remove_flip$flip)] = -gwas_pqtl_df$BETA[which(remove_flip$flip)]
    gwas_pqtl_df$EAF[which(remove_flip$flip)] = 1-gwas_pqtl_df$EAF[which(remove_flip$flip)]
    gwas_pqtl_df = gwas_pqtl_df[remove_flip$keep,]
    gwas_df = gwas_pqtl_df %>% select(SNP,CHR,A1=A1.y,A2=A2.y,freq=EAF,b=BETA,se=SE,p=P,N=N.x)
    pqtl_df = gwas_pqtl_df %>% select(SNP,CHR,A1=A1.y,A2=A2.y,freq,b,se,p,N=N.y) 

    pqtl_fn = paste0('./pqtl/',s,'.txt')
    fwrite(pqtl_df,pqtl_fn,sep='\t')
    gwas_fn = paste0('./gwas/',s,'.txt')
    fwrite(gwas_df,gwas_fn,sep='\t')
    write.table(pqtl_df[,1],paste0('./pqtl/',s,'.snp'),quote=F,row.names=F,col.names=F)
    prefix = paste0('/home/panwei/shared/UKBiobankIndiv/imputed/pgen/ukbb_chr',chr,'_1')
    plink_command = paste0("module load plink/2.00-alpha-091019; \n plink2 --pfile ", prefix, " 'vzs' --keep ~/UKB_GR_whiteUnrel.id --extract ./pqtl/",s , ".snp ", 
        " --maf 0.01 --hwe 1e-6 --geno 0.1 --make-bed --out ./LDref/", s, "_ref" )
    p = system(plink_command,intern=TRUE)

    pqtl_df = fread(paste0('./pqtl/',s,'.txt'))
    pqtl_fn = paste0('./pqtl/cojo/',s,'.txt')
    fwrite(pqtl_df[,-2],pqtl_fn,sep='\t')
    cojo_pqtl_cm = paste0("~/gcta_1.92.3beta3/gcta64 --bfile ./LDref/", s , "_ref  --chr ", chr," --cojo-file ", pqtl_fn, " --cojo-slct --cojo-p 5e-6 --out ./pqtl/cojo/", s)
    system(cojo_pqtl_cm)
    if(!file.exists(paste0('./pqtl/cojo/', s,'.jma.cojo'))){print(paste0('Skip ',s));next;}
    pqtl_cojo_res = fread(paste0('./pqtl/cojo/', s,'.jma.cojo'))
    write.table(pqtl_cojo_res$SNP,paste0('./pqtl/cojo/',s,'.cojo.snp'),quote=F,row.names=F,col.names=F)
    if(length(pqtl_cojo_res$SNP)<3){print(paste0('Skip ',s));next;}

    gwas_df = fread(paste0('./gwas/',s,'.txt'))
    gwas_fn = paste0('./gwas/cojo/',s,'.txt')
    fwrite(gwas_df %>% select(-CHR) ,gwas_fn,sep='\t')
    cojo_gwas_cm = paste0("~/gcta_1.92.3beta3/gcta64 --bfile ./LDref/", s , "_ref  --chr ", chr," --cojo-file ", gwas_fn, " --cojo-slct --cojo-p 5e-6 --out ./gwas/cojo/", s)
    system(cojo_gwas_cm)
    if(file.exists(paste0('./gwas/cojo/',s,'.jma.cojo'))){
        gwas_cojo_res = fread(paste0('./gwas/cojo/', s,'.jma.cojo'))
        IV = union(pqtl_cojo_res$SNP,gwas_cojo_res$SNP)
    }else{
        IV = pqtl_cojo_res$SNP
    }
    write.table(IV,paste0('./pqtl/',s,'.IV'),quote=F,row.names=F,col.names=F)
    IV = fread(paste0('./pqtl/',s,'.IV'),header=F)$V1
    joint_pqtl_cm = paste0("~/gcta_1.92.3beta3/gcta64 --bfile ./LDref/", s , "_ref  --chr ", chr," --cojo-file ", pqtl_fn, " --extract ./pqtl/", s,".IV --cojo-joint  --out ./pqtl/cojo/", s)
    p = system(joint_pqtl_cm,intern=F)
    if(p==1){print(paste0('Collinear-',s)); next;}

    LD_mat = fread(paste0('./pqtl/cojo/', s,'.ldr.cojo'))
    LD_mat = LD_mat[,2:(ncol(LD_mat)-1)]
    LD_mat = as.matrix(LD_mat)
    rownames(LD_mat) = colnames(LD_mat)

    pqtl_df = pqtl_df %>% filter(SNP %in% IV)
    pqtl_df = pqtl_df[match(colnames(LD_mat),pqtl_df$SNP),]
    pqtl_df$p = as.numeric(pqtl_df$p)
    pqtl_df$z = pqtl_df$b/pqtl_df$se
    pqtl_cor = pqtl_df$b / sqrt(pqtl_df$b^2 + (pqtl_df$N - 2) * pqtl_df$se^2)
    pqtl_df$cor = pqtl_cor
    pqtl_df$se_cor = pqtl_df$se * pqtl_df$cor/pqtl_df$b
    pqtl_df$bJ = solve(LD_mat) %*% pqtl_cor
    # sigma2_J = 1 - 2 * pqtl_cor %*% pqtl_df$bJ + t(pqtl_df$bJ) %*% LD_mat %*% pqtl_df$bJ
    # bJ_pqtl_var = solve(LD_mat)/median(pqtl_df$N) * as.numeric(sigma2_J)
    # bJ_pqtl_var1 = solve(LD_mat) %*% (pqtl_df$se_cor %o% pqtl_df$se_cor * LD_mat) %*% solve(LD_mat)

    gwas_df = gwas_df %>% filter(SNP %in% IV)
    gwas_df = gwas_df[match(colnames(LD_mat),gwas_df$SNP),]
    gwas_df$p = as.numeric(gwas_df$p)
    gwas_df$z = gwas_df$b/gwas_df$se
    gwas_cor = gwas_df$b / sqrt(gwas_df$b^2 + (gwas_df$N - 2) * gwas_df$se^2)
    gwas_df$cor = gwas_cor
    gwas_df$se_cor = gwas_df$se * gwas_df$cor/gwas_df$b
    gwas_df$bJ = solve(LD_mat) %*% gwas_cor
    # sigma2_J = 1 - 2 * gwas_cor %*% gwas_df$bJ + t(gwas_df$bJ) %*% LD_mat %*% gwas_df$bJ
    # bJ_gwas_var = solve(LD_mat)/median(gwas_df$N) * as.numeric(sigma2_J)
    # bJ_gwas_var1 = solve(LD_mat) %*% (gwas_df$se_cor %o% gwas_df$se_cor * LD_mat) %*% solve(LD_mat)

    mr_dat = list(b_exp=pqtl_df$bJ, b_out=gwas_df$bJ, N1=median(pqtl_df$N), N2=median(gwas_df$N), LD_mat = LD_mat,
        exp_df = pqtl_df, out_df = gwas_df,
        # Sig_exp1 = bJ_pqtl_var1, Sig_out1 = bJ_gwas_var1,
        exp_IV=pqtl_cojo_res$SNP,out_IV=setdiff(IV,pqtl_cojo_res$SNP))
    save(mr_dat, file=paste0('./MRdat/',s,'.RData'))


}
