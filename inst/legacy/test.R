simGO<-function(prefix_bin=NULL, #currently only consider .bed .bim .bam format data 
                #for different ancestries, prefix_bin is a vector 
                
                h2=0.4,       #heratability of trait, accounted by genotype dataset and omics dataset 
                h2snps=0.2,  #heratability accounted by genotype data 
                h2omics=0.7, #heratability of omics data,  s2omics_g/s2omic					  
                s2beta_omics=0.01,     #variance for the beta_omics effect 
                s2m=2, # variance of each omics measurement
                
                #for simplify, we assume all the snps 
                nSNPs_causal=NULL,   #3000, number of causal QTLs
                nSNPs_causal_per=NULL, #1%, percentage of causal SNPs
                nOmics_features=1000, # number of omics features 					  
                nSNPs_pleiotropy=1000, # number of causal SNPs shared 
                #nSNPs_pleiotropy_per=NULL, # percentage of causal SNPs shared 
                
                snp_cor=0.5, #genetic correlation based on genotype data 
                gene_cor=NULL, #correlation of gene expression data within a specified region
                #currently, we only considered the gene correlation within a specified region 
                
                #Parameters for SNPs
                SNPs_archi="Sparse", #Infinite,Sparse,BSLMM,BayesR
                SNPs_var=NULL,#Variances of SNPs within each group,
                #c(1),Infinite; c(1,0), Sparse; c(0.8,0.2), BSLMM assumption; c(1,0.001,0.01,0) , BayesR assumption
                # SNPs with non-zero variances are regarded as causal SNPs              
                pSNPs_var=NULL,#proportions of SNPs with specified variance within each group,correspond to SNPs_var, 
                #c(0.2,0.8), BSLMM assumption; c(0.05,0.05,0.1,,0.8), BayesR assumption 
                nSNPs_var=NULL,#number of SNPs with specified variance within each group,correspond to SNPs_var, 
                SNPs_causal_mode="Random", #Random, means random sampling, ignore LD
                #Random_uniform, means random sampling within the even blocks, ignore LD
                #Random_blocks, means random sampling within blocks,make sure at least 1 SNP within each block
                #Exact_blocks, means extract all SNPs within the specified blocks, can be used to extract cis-SNPs within gene
                
                SNPs_causal_blocks=NULL, #sampling SNPs from the provided LD blocks
                #dataframe, 1st is blocks start point, 2st is blocks end point
                sSNPs_causal=1, #{√(2pq)}^s, control the degree of negative selection,
                #genotype matrix should not be scaled in corresponding to selection
                # 0 means no negative selection, 1 means strong negative , 
                # 0.25 means mild negative selection, -1 means positive selection
                #genetic correlation based on SNPs data
                gen_cor=NULL, #genetic correlation based on genotype data, off-diagonal value of genetic correlation matrix 
                nSNPs_causal_shared=NULL, # number of shared causal SNPs across all causal SNPs
                pSNPs_causal_shared=0.2, #percentage of shared causal SNPs across all causal SNPs
                
                #Omics parameters 
                nOmics_features=1000, # number of omics features 		
                nOmics_causal=500, #number of omics features contribute to phenotype
                pOmics_causal=0.5, #percentnage of omics features contribute to phenotype
                
                Omics_SNPs_archi="Infinite", #Infinite,Sparse,BSLMM,BayesR    
                #within each region, how SNPs contributes to omics(e.g., gene expression)   
                #Infinite, all causal SNPs contribute to omics
                #BSLMM, small proportions of SNPs with large effects, and large proportions of SNPs with small effects contributes to omics features
                #BayesR, similar to BSLMM but with more groups
                
                pOmics_SNPs_causal=0.1, #percentage of causal SNPs contributes to omics features
                nOmics_SNPs_causal=NULL, #number of causal SNPs contributes to omics features
                Omics_SNPs_causal_mode="Random", #how to select causal SNPs whereas contributing to Omics features,
                #Random, means random sampling,  set of causal SNPs contributes to each Omics feature will be different
                #Random_blocks, means random sampling within blocks, make sure at least 1 SNP within each block
                #Exact_blocks, means extract all SNPs within the specified blocks, can be used for gene expression
                Omics_SNPs_causal_blocks=NULL, #sampling SNPs from the provided LD blocks
                #
                
                ##correlated omics features (e.g., gene expression)
                
                
                
                
                
                write_data=TRUE, #write the generated phe and omics data as file 
                output_path=getwd() # default is the current path
){
  
  #check the input of causal SNPs 
  if(is.null(nSNPs_causal)&is.null(nSNPs_causal_per))stop("nSNPs_causal and nSNPs_causal_per could be null simultaneously!")
  if(!is.null(nSNPs_causal)&!is.null(nSNPs_causal_per))stop("nSNPs_causal and nSNPs_causal_per could be specified simultaneously!")
  if(!is.null(nSNPs_causal_per)){
    nSNPs_causal=ceiling(nSNP*nSNPs_causal_per,0);cat(paste0(nSNPs_causal_per*100,"% SNPs are setted as casusal SNPs in the simulation! \n"))
  }else{
    cat(paste0(nSNPs_causal," SNPs are setted as casusal SNPs in the simulation! \n"))
  }
  
  #SNPs genetic architectures 
  if(SNPs_archi=="Infinite"){
    
  }else if(SNPs_archi=="Sparse"){
    
  }else if(SNPs_archi=="BSLMM"){
    
  }else if(SNPs_archi=="BayesR"){
    
  }
  
  
  
  
  #generate the genetic correlation 
  if(!is.null(gen_cor)){
    beta=create_gen_cor(nSNPs_causal,nSNPs_causal_shared,snp_cor,seed)
  }
  
  
  
  #define the heratability and the variance components 		 
  # s2y=s2omics+s2snps+s2e #sigma2_y=sigma2_omics+sigma2_snps+sigma2_e 
  
  # s2omics=s2omics_g+s2omics_ng #sigma2_omics=sigma2_omics_genetics+sigma2_snps+sigma2_omics_non_genetics 
  
  # h2snps=s2snps/s2y
  # h2omics=s2omics_g/s2omics #omics effects accounted by genetics 
  
  # rs2omics=s2omics/s2y   #ratio of omics to total phenotype variance 
  
  # h2=h2omics*rs2omics+h2snps=(s2omics_g/s2omics)*(s2omics/s2y)+s2snps/s2y=(s2omics_g+s2snps)/s2y #heritability of trait 
  
  #initialize of differnt parameters 
  # h2=0.4 
  # h2snps=0.2 
  # h2omics=0.7
  # rs2omics=(h2-h2snps)/h2omics
  # s2m=2 # variance of each omics measurement
  # nOmics_features=1000 # number of omics features 
  # s2beta_omics=0.01     #variance for the beta_omics effect 
  
  ##for simplify, we assume all the snps 
  # nCausal=3000; #10% causal
  # pSNPs_causal=1     #percentage of causal contributes to snps 
  # pOmics_causal=0.8  #percentage of causal contributes to omics features 
  ##				   here we assume the percentage is the same for all omics features 
  
  # equal_Omics_causal=FALSE # if the same set of causals contributes to each Omics data
  
  
  setMKLthreads(4)
  library(data.table)
  #write data path
  if(!file.exists(output_path)){
    
    dir.create(output_path, recursive = TRUE)
    
  }
  
  rs2omics=(h2-h2snps)/h2omics # the ratio of omics variance to total phenotypic variance 
  cat(paste0("The ratio of omics data to total phenotypic variance is ",round(rs2omics,3),"!\n"))
  
  
  
  ############read genotype data 
  #Hapmap3 
  snp = bed_to_numeric(prefix_bin)
  bim = read.table(paste0(prefix_bin,".bim"),header=F)
  fam = read.table(paste0(prefix_bin,".fam"),header=F)
  
  rownames(snp) = as.character(fam[, 1])
  colnames(snp) = as.character(bim[, 2])
  nSNP=ncol(snp);
  nIND=nrow(snp);
  snp=scaleM_Y(snp)  #scale and center snp data
  ############
  
  #constant heritability for each SNP
  # Sample the 3000 variants
  set.seed(sim_seed)
  nOmics_causal=nSNPs_causal*pOmics_causal
  pos_causal=sample(seq(1, nSNP), nSNPs_causal)
  pos_SNPs_causal = sample(pos_causal, nSNPs_causal)
  
  snp_set=snp[,pos_SNPs_causal] #scaleM_Y for scaling and centring genotype matrix
  
  if(nOmics_causal==nSNPs_causal&equal_Omics_causal){ # all causal SNPs contributes to the omics data 
    
    pos_Omics_causal=pos_Omics_causal
    omics_snp_set=snp_set #scaleM_Y for scaling and centring genotype matrix
    
  }else{
    
    if(equal_Omics_causal){ # the same set of causals contributes to each Omics data
      
      pos_Omics_causal   = sample(pos_causal, nOmics_causal)
      omics_snp_set=snp[,pos_Omics_causal]#/sqrt(ncol(nOmics_causal))
      
    }else{
      
      pos_Omics_causal   =matrix(NA,ncol=nOmics_features,nrow=nOmics_causal)
      for (i in 1:nOmics_features){
        pos_Omics_causal[,i]=sample(pos_causal, nOmics_causal)
      }
    }	
    
    
  }
  
  
  #sample omics effects, snps effects and omics weights
  weight_omics =matrix(NA,ncol=nOmics_features,nrow=nOmics_causal)
  set.seed(sim_seed+1) 
  beta_omics=rnorm(nOmics_features,0,sqrt(s2beta_omics))
  beta_snps=rnorm(nSNPs_causal)
  
  for (i in 1:nOmics_features){
    weight_omics[,i]=rnorm(nOmics_causal)  
  }
  
  
  # simulate omics measurement for each features 
  # omics_measurement = yomics_g + yomics_ng
  yomics_g=matrix(NA,ncol=nOmics_features,nrow=nIND) #measurement of omics data which accounted by genetics
  yomics_ng=matrix(NA,ncol=nOmics_features,nrow=nIND) #measurement of omics data which accounted by non-genetics
  cat(paste0("Start simulating ",nOmics_features, " Omics-features!\n"))
  for (i in 1:nOmics_features){ #for each omics features 
    
    cat(paste0("Simulated Omics-features ",i," done!\n"))
    if(!equal_Omics_causal){
      
      omics_snp_set=snp[,pos_Omics_causal[,i]]
    }
    
    yomics_g[,i]=omics_snp_set%*%weight_omics[,i]/sqrt(nOmics_causal) #measurement of omics data which accounted by genetics
    yomics_g[,i]=sqrt(s2m*h2omics)*yomics_g[,i]#
    varG=var(yomics_g[,i])
    varE=s2m*(1-h2omics)
    yomics_ng[,i]=rnorm(nIND,0,sd=sqrt(varE))
    varE=var(yomics_ng[,i])
    #cat("nOmics_features=",i,"varG= ",varG,"varE= ",varE,"varM= ",varG+varE,"h2omics ",varG/(varG+varE),"\n")
  }
  
  #M matrix, the measurement of each omics features
  yomics=yomics_g+yomics_ng
  
  uomics_g = as.vector(yomics_g%*%beta_omics) #genetic part contributed by omics, individual level
  uomics_ng   = as.vector(yomics_ng%*%beta_omics)    #residual part contributed by omics , individual level
  uomics = uomics_g + uomics_ng   # omics effect estiamted by REML, individual level
  
  # cat("computing required numbers\n")
  s2omics=var(uomics)
  # cat("s2omics= ",s2omics,"\n")
  s2y=h2omics*s2omics/(h2-h2snps)      # s2omics/rs2omics, s2y: total phenotypic variance
  # cat("s2y= ",s2y,"\n")
  s2snps=h2snps*s2y
  # cat("s2snps= ",s2snps,"\n")
  rs2omics=(h2-h2snps)/h2omics
  # cat("rs2omics= ",rs2omics,"\n")
  s2e=(1-rs2omics-h2snps)*s2y
  # cat("s2e= ",s2e,"\n")
  
  # genotype of SNPs 
  usnps=as.vector(snp_set%*%beta_snps)
  usnps=sqrt(s2snps)*usnps/sd(usnps)
  epsilon=rnorm(nIND,0,sqrt(s2e))
  y=uomics+usnps+epsilon
  
  gcta_phe=cbind(fam[,1:2],y)
  real_phe=cbind(fam[,1:2],usnps+uomics_g)
  
  #write data 
  setwd(output_path)
  if(write_data){
    
    cat(paste0("Simulated data are saved in path:",getwd(),"!\n"))
    fwrite(data.table(gcta_phe),"sim_pheno.txt",row.names=F,col.names=F,quote=F,sep=" ")
    fwrite(data.table(yomics),"sim_OmicsM.txt", row.names=F,col.names=F,quote=F,sep=" ")	
    fwrite(data.table(real_phe),"sim_tgv.txt",row.names=F,col.names=F,quote=F,sep=" ")
    tmp=yomics%*%t(yomics)+diag(0.0001,nrow(yomics)); #make sure tmp is invertible 
    Gm=lotri_matrix_col3(tmp/mean(diag(tmp)),nOmics_causal);
    rm(tmp,yomics);gc();
    fwrite(data.table(Gm),"sim_Omics.grm", row.names=F,col.names=F,quote=F,sep=" ");	
    system(paste0("gzip sim_Omics.grm"))
    fwrite(data.table(fam[,1:2]),"sim_Omics.grm.id", row.names=F,col.names=F,quote=F,sep=" ");			
  }	
  
}
