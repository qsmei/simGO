#generate genetic correaltion for multi traits, multi ancestries, 
#and multi omics features (e.g., correlation between gene expression)
create_gen_cor<-fucntion(
                         nSNPs_causal=c(2000,2000), #number of causal SNPs, 
                         nSNPs_causal_shared=1000, # number of shared causal SNPs
                         nSNPs_causal_shared_per=0.5, #percentage of shared causal SNPs
                         h2_SNPs=c(0.5,0.5), #snps heritability
                         gen_cor=0.8, #genetic correlations for different traits ()
                         seed=1998 #heritability for each trait
                         ){
  
  set.seed(seed)
  #for two traits situation
  nSNPs_causal_non_shared=nSNPs_causal-nSNPs_causal_shared #for different tratis,  the non share snps could be different
  gen_cor_mat=diag(h2_SNPs/nSNPs_causal)
  gen_cor_mat[lower.tri(gen_cor_mat)]=gen_cor*sqrt(prod(h2_SNPs))/nSNPs_causal_shared
  #generated shared beta effects
  ####beta_snp_shared <- MASS::mvrnorm(nSNPs_causal_shared,c(0,0),gen_cor_mat) #assumed mean effects for traits 
  #### Cholesky decomposition of the covariance matrix
  chol_cov <- chol(gen_cor_mat)
  # Generate n samples of standard normal random variables (p variables)
  z <- matrix(rnorm(nSNPs_causal_shared * 2), nrow = nSNPs_causal_shared, ncol = 2)
  beta_snp_shared=(matrix(c(0,0), nrow = nSNPs_causal_shared, ncol = 2, byrow = TRUE) + z %*% chol_cov)
  #####
  
  #non shared beta coefficients
  beta_nonshared_1 <- rnorm(nSNPs_causal_non_shared[1],0,sd=sqrt(1-h2[1]))
  beta_nonshared_2 <- rnorm(nSNPs_causal_non_shared[2],0,sd=sqrt(1-h2[2]))
  
  beta1=rbind(beta_snp_shared[,1],beta_nonshared_1)
  beta2=rbind(beta_snp_shared[,2],beta_nonshared_2)
  
}
  
  
  
  
  
  
}