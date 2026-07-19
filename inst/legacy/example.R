#same ancestry,Sparse
nSNPs=10000
nTraits=2
SNPs_archi="Sparse" #InfiniteSparseBSLMMBayesR
SNPs_causal_var=1#Variances of SNPs within each group
pSNPs_causal_var=0.1#proportions of SNPs with specified variance within each groupcorrespond to SNPs_var 
nSNPs_causal_var=NULL#number of SNPs with specified variance within each groupcorrespond to SNPs_var 
SNPs_causal_mode="Random" #Random means random sampling ignore LD
SNPs_causal_blocks=NULL #sampling SNPs from the provided LD blocks
sSNPs_causal=1 #{√(2pq)}^s control the degree of negative selection
SNPs_cor=0.8 #genetic correlation based on genotype data off-diagonal value of genetic correlation matrix 
nSNPs_causal_shared=400 # number of shared causal SNPs across all causal SNPs
#pSNPs_causal_shared=0.2 #percentage of shared causal SNPs across all causal SNPs

        simGO(prefix_bin=NULL,nTraits=2,
                h2SNPs=list(0.4,0.8),  #heratability accounted by genotype data 
                SNPs_archi="Sparse", #Infinite,Sparse,BSLMM,BayesR
                SNPs_causal_var=1,#Variances of SNPs within each group,
                pSNPs_causal_var=NULL,#proportions of SNPs with specified variance within each group,correspond to SNPs_var, 
                nSNPs_causal_var=1000,#number of SNPs with specified variance within each group,correspond to SNPs_var, 
                #if input of SNPs(pSNPs_causal_var <1 or nSNPs_causal_var < nSNPs),the rest SNPs are stand for zero-effect SNPs 
                SNPs_causal_mode="Random", #Random, means random sampling, ignore LD
                sSNPs_negative=-1, #{√(2pq)}^s, control the degree of negative selection,
                SNPs_cor=NULL, #genetic correlation based on genotype data, off-diagonal value of genetic correlation matrix
                nSNPs_causal_shared=1000, # number of shared causal SNPs across all causal SNPs
                Res_cor=0.2,
                write_data=TRUE, #write the generated phe and omics data as file 
                output_path=getwd() # default is the current path
          )


#same ancestry,Sparse
nSNPs=10000
nTraits=2
h2=0.4
h2SNPs=list(0.4,0.8)
h2Omics=NULL
SNPs_archi="Sparse" #InfiniteSparseBSLMMBayesR
SNPs_causal_var=1#Variances of SNPs within each group
#pSNPs_causal_var=c(0.3,0.7)#proportions of SNPs with specified variance within each group correspond to SNPs_var 
#means 3000 SNPs with large effects, 7000 SNPs with small effects
nSNPs_causal_var=1000#number of SNPs with specified variance within each groupcorrespond to SNPs_var 
SNPs_causal_mode="Random" #Random means random sampling ignore LD
SNPs_causal_blocks=NULL #sampling SNPs from the provided LD blocks
sSNPs_causal=1 #{√(2pq)}^s control the degree of negative selection
SNPs_cor=0.8 #genetic correlation based on genotype data off-diagonal value of genetic correlation matrix 
nSNPs_causal_shared=1000 # number of shared causal SNPs across all causal SNPs
# For trait1,
# 1000 SNPs within 3000 large effects SNPs are shared across traits 
# 2000 SNPs within 7000 small effects SNPs are shared across traits 
# For trait2,
# 1000 SNPs within 2000 large effects SNPs are shared across traits 
# 2000 SNPs within 8000 small effects SNPs are shared across traits 
pSNPs_causal_shared=NULL #percentage of shared causal SNPs across all causal SNPs
nSNPs=list(10000,10000)
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
  sSNPs_causal=check_input_multi(nTraits,sSNPs_causal)
  nSNPs=check_input_multi(nTraits,nSNPs) #input SNPs 
}




#same ancestry,BSLMM
nSNPs=10000
nTraits=2
h2=0.4
h2SNPs=list(0.4,0.8)
h2Omics=NULL
SNPs_archi="BSLMM" #InfiniteSparseBSLMMBayesR
SNPs_causal_var=c(0.8,0.2)#Variances of SNPs within each group
#pSNPs_causal_var=c(0.3,0.7)#proportions of SNPs with specified variance within each group correspond to SNPs_var 
                           #means 3000 SNPs with large effects, 7000 SNPs with small effects
nSNPs_causal_var=list(c(3000,7000),c(3000,7000))#number of SNPs with specified variance within each groupcorrespond to SNPs_var 
SNPs_causal_mode="Random" #Random means random sampling ignore LD
SNPs_causal_blocks=NULL #sampling SNPs from the provided LD blocks
sSNPs_causal=1 #{√(2pq)}^s control the degree of negative selection
SNPs_cor=0.8 #genetic correlation based on genotype data off-diagonal value of genetic correlation matrix 
nSNPs_causal_shared=c(2000,6000) # number of shared causal SNPs across all causal SNPs
                                 # For trait1,
                                 # 1000 SNPs within 3000 large effects SNPs are shared across traits 
                                 # 2000 SNPs within 7000 small effects SNPs are shared across traits 
                                 # For trait2,
                                 # 1000 SNPs within 2000 large effects SNPs are shared across traits 
                                 # 2000 SNPs within 8000 small effects SNPs are shared across traits 
pSNPs_causal_shared=NULL #percentage of shared causal SNPs across all causal SNPs
nSNPs=list(10000,10000)
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
  sSNPs_causal=check_input_multi(nTraits,sSNPs_causal)
  nSNPs=check_input_multi(nTraits,nSNPs) #input SNPs 
}


#same ancestry,BayesR
nSNPs=10000
nTraits=2
h2=0.4
h2SNPs=list(0.4,0.8)
h2Omics=NULL
SNPs_archi="BayesR" #InfiniteSparseBSLMMBayesR
SNPs_causal_var=c(1,0.01,0.001)#Variances of SNPs within each group
#pSNPs_causal_var=c(0.3,0.7)#proportions of SNPs with specified variance within each group correspond to SNPs_var 
#means 3000 SNPs with large effects, 7000 SNPs with small effects
nSNPs_causal_var=c(1000,2000,6000)#number of SNPs with specified variance within each groupcorrespond to SNPs_var 
#nSNPs_causal_var=list(c(3000,7000),c(4000,8000))#number of SNPs with specified variance within each groupcorrespond to SNPs_var 
SNPs_causal_mode="Random" #Random means random sampling ignore LD
SNPs_causal_blocks=NULL #sampling SNPs from the provided LD blocks
sSNPs_causal=1 #{√(2pq)}^s control the degree of negative selection
SNPs_cor=0.8 #genetic correlation based on genotype data off-diagonal value of genetic correlation matrix 
nSNPs_causal_shared=c(1000,2000,6000) # number of shared causal SNPs across all causal SNPs
# For trait1,
# 1000 SNPs within 3000 large effects SNPs are shared across traits 
# 2000 SNPs within 7000 small effects SNPs are shared across traits 
# For trait2,
# 1000 SNPs within 2000 large effects SNPs are shared across traits 
# 2000 SNPs within 8000 small effects SNPs are shared across traits 
pSNPs_causal_shared=NULL #percentage of shared causal SNPs across all causal SNPs
nSNPs=list(10000,10000)
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
  sSNPs_causal=check_input_multi(nTraits,sSNPs_causal)
  nSNPs=check_input_multi(nTraits,nSNPs) #input SNPs 
}
