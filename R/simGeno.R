setwd("/restricted/projectnb/fhs_data/qsmei/Phase3_1000G_imputed/1000GP_Phase3/Hapmap3_1000G/simulation/EUR/Rcpp_code")
Rcpp::sourceCpp("readbin.cpp")
Rcpp::sourceCpp("scaleM.cpp")
simGO<-function(prefix_bin=NULL, #currently only consider .bed .bim .bam format data 
                #for different ancestries, prefix_bin is a vector 
                nTraits=1,    #number of simulated traits
                h2=0.4,       #heratability of trait, accounted by genotype dataset and omics dataset 
                h2SNPs=0.2,  #heratability accounted by genotype data 
                h2Omics=0.7, #heratability of omics data,  s2omics_g/s2omic					  
                #s2beta_omics=0.01,     #variance for the beta_omics effect 
                #s2m=2, # variance of each omics measurement
                
                snp_cor=0.5, #genetic correlation based on genotype data 
                gene_cor=NULL, #correlation of gene expression data within a specified region
                #currently, we only considered the gene correlation within a specified region 
                
                #Parameters for SNPs
                #for multiple traits scenario, if class of (SNPs_archi,SNPs_causal_var,pSNPs_causal_var....)!=list, means all traits shared the same structure
                SNPs_archi="Sparse", #Infinite,Sparse,BSLMM,BayesR
                SNPs_causal_var=1,#Variances of SNPs within each group,
                #c(1),Infinite; c(1,0), Sparse; c(0.8,0.2), BSLMM assumption; c(1,0.001,0.01,0) , BayesR assumption
                # SNPs with non-zero variances are regarded as causal SNPs              
                pSNPs_causal_var=NULL,#proportions of SNPs with specified variance within each group,correspond to SNPs_var, 
                #c(0.2,0.8), BSLMM assumption; c(0.05,0.05,0.1,,0.8), BayesR assumption
                
                nSNPs_causal_var=NULL,#number of SNPs with specified variance within each group,correspond to SNPs_var, 
                #if input of SNPs(pSNPs_causal_var <1 or nSNPs_causal_var < nSNPs),the rest SNPs are stand for zero-effect SNPs 
                SNPs_causal_mode="Random", #Random, means random sampling, ignore LD
                #Random_uniform, means random sampling within the even blocks, ignore LD
                #Random_blocks, means random sampling within blocks,make sure at least 1 SNP within each block
                #Exact_blocks, means extract all SNPs within the specified blocks, can be used to extract cis-SNPs within gene
                
                SNPs_causal_blocks=NULL, #sampling SNPs from the provided LD blocks
                #dataframe, 1st is blocks start point, 2st is blocks end point
                sSNPs_negative=1, #{√(2pq)}^s, control the degree of negative selection,
                #genotype matrix should not be scaled in corresponding to selection
                # 0 means no negative selection, 1 means strong negative , 
                # 0.25 means mild negative selection, -1 means positive selection
                #genetic correlation based on SNPs data
                SNPs_cor=NULL, #genetic correlation based on genotype data, off-diagonal value of genetic correlation matrix
                pSNPs_causal_same_sets=1, # 1 means different traits have same causal SNPs set, under infinite situation, all traits have the same causal SNPs set
                # 0 means different traits have different causal SNPs set
                
                nSNPs_causal_shared=NULL, # number of shared causal SNPs across all causal SNPs
                pSNPs_causal_shared=0.2, #percentage of shared causal SNPs across all causal SNPs
                #for nSNPs_causal_shared and pSNPs_causal_shared, they should be a vector correspond to SNPs_causal_var for BSLMM and BayesR,
                #e.g. In total 2000 shared SNPs, 500 from large effects, 1500 from small effects
                
                Res_cor=0.2,
                
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
                
                seed=1998,
                write_data=TRUE, #write the generated phe and omics data as file 
                output_path=getwd() # default is the current path
){
  #setMKLthreads(4)
  library(data.table)
  #write data path
  if(!file.exists(output_path)){
    dir.create(output_path, recursive = TRUE)
  }
  ############read genotype data 
  #.bed format data 
  SNPs_mat = bed_to_numeric(prefix_bin)
  bim = read.table(paste0(prefix_bin,".bim"),header=F)
  fam = read.table(paste0(prefix_bin,".fam"),header=F)
  rownames(SNPs_mat) = as.character(fam[, 1])
  colnames(SNPs_mat) = as.character(bim[, 2])
  nSNPs=ncol(SNPs_mat);
  nINDs=nrow(SNPs_mat);
  
  p=colMeans(SNPs_mat)/2
  sr_2pq=(2*p*(1-p))^{0.5*sSNPs_negative} #square root of 2pq, s=0 means no negetive selection, s=-1 means strong negative selection, s= -0.5 for mild negative selection
  #snp=scaleM_Y(snp)  #scale and center snp data
  #consider negative selection as BayesS, thus not scale genotype
  
  ############
  #### SNPs
  ############
  if(is.null(nSNPs_causal_var)&is.null(pSNPs_causal_var))stop("nSNPs_causal_var and pSNPs_causal_var could not be null simultaneously!")
  if(!is.null(nSNPs_causal_var)&!is.null(pSNPs_causal_var))stop("nSNPs_causal_var and pSNPs_causal_var could not be specified simultaneously!")
  
  #Inputs of multiple traits could be list or vector, if nTraits>1, and input parameter is not vector, this means all traits shared the same parameters
  if(nTraits>=2){
    h2=check_input_multi(nTraits,h2) 
    h2SNPs=check_input_multi(nTraits,h2SNPs) 
    h2Omics=check_input_multi(nTraits,h2Omics) 
    SNPs_archi=check_input_multi(nTraits,SNPs_archi) 
    SNPs_causal_var=check_input_multi(nTraits,SNPs_causal_var) 
    nSNPs_causal_var=check_input_multi(nTraits,nSNPs_causal_var) 
    pSNPs_causal_var=check_input_multi(nTraits,pSNPs_causal_var) 
    SNPs_causal_mode=check_input_multi(nTraits,SNPs_causal_mode) 
    SNPs_causal_blocks=check_input_multi(nTraits,SNPs_causal_blocks) 
    sSNPs_negative=check_input_multi(nTraits,sSNPs_negative)
    nSNPs=check_input_multi(nTraits,nSNPs) #input SNPs 
    nINDs=check_input_multi(nTraits,nINDs)
  }
  #check the input of causal SNPs 
  nSNPs_causal_var=check_input_SNPs(nTraits,SNPs_archi,nSNPs,pSNPs_causal_var,nSNPs_causal_var)
  
  cat(paste0("Sampling beta effects for SNPs....\n"))
  #sample beta effects for SNPs across traits
  beta_SNPs=beta_sample_SNPs(nTraits,SNPs_archi,SNPs_causal_mode,nSNPs,h2SNPs,nSNPs_causal_var,SNPs_cor,nSNPs_causal_shared,seed)
  cat(paste0("Sampling beta effects for SNPs done!\n"))
  #Assign sampled beta effects into SNPs
  
  
  #sample error terms for Individuals across traits
  cat(paste0("Sampling residual effects for INDs....\n"))
  residual_INDs=sample_residual(nTraits,Res_cor,h2SNPs,nINDs,seed)
  cat(paste0("Sampling residual effects for INDs done!\n"))
  
  #create phenotype 
  cat(paste0("Generating phenotype for INDs....\n"))
  sim_phe=matrix(NA,nrow=nINDs[[1]],ncol=nTraits+2)
  sim_g=matrix(NA,nrow=nINDs[[1]],ncol=nTraits+2)
  sim_phe[,1]=sim_g[,1]=fam[,1];sim_phe[,2]=sim_g[,2]=fam[,2]
  
  colnames(sim_phe)=c("Fam","IND",paste0("Phe_Trait",1:nTraits))
  colnames(sim_g)=c("Fam","IND",paste0("Genetic_Trait",1:nTraits))
  
  for(i in 1:nTraits){
    
    # genotype of SNPs 
    i_index=as.numeric(beta_SNPs[,i*3-2])
    
    usnps=as.vector(SNPs_mat[,i_index]%*%as.numeric(beta_SNPs[,i*3]))
    sim_g[,i+2]=usnps
    sim_phe[,i+2]=usnps+residual_INDs[,i]
  }
  cat(paste0("Generating phenotype for INDs done!\n"))
  
  #write data 
  setwd(output_path)
  if(write_data){
    
    cat(paste0("Simulated data are saved in path:",getwd(),"!\n"))
    fwrite(data.table(sim_phe),"sim_phe.txt",row.names=F,col.names=T,quote=F,sep=" ")
    fwrite(data.table(sim_phe),"sim_g.txt",row.names=F,col.names=T,quote=F,sep=" ")
    fwrite(data.table(beta_SNPs),"sim_beta.txt",row.names=F,col.names=T,quote=F,sep=" ")
  }	
  return(list(sim_phe=sim_phe,sim_g=sim_g,sim_beta=beta_SNPs))
}
