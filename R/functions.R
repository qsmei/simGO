#convert single input as multiple input under multiple traits scenario
check_input_multi<-function(nTraits=2,input=NULL){
  if(class(input)=="list"&length(input)!=nTraits)stop(paste0("The length of input para:",deparse(substitute(input))," doesn't match the length of Traits!"))
  if(class(input)!="list"){input=replicate(nTraits, input, simplify = FALSE)};
  return(input)
}

#check the input of causal SNPs parameters
check_input_SNPs<-function(nTraits=1,SNPs_archi=NULL,nSNPs=NULL,pSNPs_causal_var=NULL,nSNPs_causal_var=NULL){
  
  nSNPs_causal=replicate(nTraits, NULL, simplify = FALSE)
  
  for(i in 1:nTraits){

    if( (!is.null(pSNPs_causal_var[[i]]) & sum(pSNPs_causal_var[[i]])!=1 ) | (!is.null(nSNPs_causal_var[[i]]) & sum(nSNPs_causal_var[[i]])!=nSNPs[[i]])){
      if(SNPs_archi[[i]]=="Infinite"|SNPs_archi[[i]]=="BSLMM"){  #Infinite means all SNPs are causal
        stop("SNPs_archi=Infinite or BSLMM  requires all SNPs are causal, please check the input:pSNPs_causal_var or nSNPs_causal_var!")
      }
    }
  
    if(SNPs_archi[[i]]=="Infinite"){
      if(!is.null(pSNPs_causal_var[[i]])){
        if(length(pSNPs_causal_var[[i]])!=1)stop("SNPs_archi=Infinite requires the length of pSNPs_causal_var equals to 1 for Trait",i,"! \n")
      }else{
        if(length(nSNPs_causal_var[[i]])!=1)stop("SNPs_archi=Infinite requires the length of nSNPs_causal_var equals to 1 for Trait",i,"! \n")
      }
      nSNPs_causal[[i]]=nSNPs[[i]]
      
    }else if(SNPs_archi[[i]]=="BSLMM"){
      if(!is.null(pSNPs_causal_var[[i]])){
        if(length(pSNPs_causal_var[[i]])!=2)stop("SNPs_archi=BSLMM requires the length of pSNPs_causal_var equals to 2 for Trait",i,"! \n")
        nSNPs_causal[[i]]=ceiling(nSNPs[[i]]*pSNPs_causal_var[[i]]);
      }else{
        if(length(nSNPs_causal_var[[i]])!=2)stop("SNPs_archi=BSLMM requires the length of nSNPs_causal_var equals to 2 for Trait",i,"! \n")
        nSNPs_causal[[i]]=nSNPs_causal_var[[i]]
      }        
      
    }else if(SNPs_archi[[i]]=="Sparse"){
      if(!is.null(pSNPs_causal_var[[i]])){
        if(length(pSNPs_causal_var[[i]])!=1)stop("SNPs_archi=Sparse requires the length of pSNPs_causal_var equals to 1 for Trait",i,"! \n")
        if(sum(pSNPs_causal_var[[i]])!=1)cat(paste0(100*(1-sum(pSNPs_causal_var[[i]]))," % SNPs are assigned zero effects for Trait",i,"! \n"))
        nSNPs_causal[[i]]=ceiling(nSNPs[[i]]*pSNPs_causal_var[[i]]);
      }else{
        if(length(nSNPs_causal_var[[i]])!=1)stop("SNPs_archi=Sparse requires the length of nSNPs_causal_var equals to 1 for Trait",i,"! \n")
        if(sum(nSNPs_causal_var[[i]])!=nSNPs[[i]])cat(paste0((nSNPs[[i]]-sum(nSNPs_causal_var[[i]]))," SNPs are assigned zero effects for Trait",i,"! \n"))
        nSNPs_causal[[i]]=nSNPs_causal_var[[i]]
      }
    }else if(SNPs_archi[[i]]=="BayesR"){
      if(!is.null(pSNPs_causal_var[[i]])){
        if(sum(pSNPs_causal_var[[i]])!=1)cat(paste0(100*(1-sum(pSNPs_causal_var[[i]]))," % SNPs are assigned zero effects for Trait",i,"! \n"))
        nSNPs_causal[[i]]=ceiling(nSNPs[[i]]*pSNPs_causal_var[[i]])
      }else{
        if(sum(nSNPs_causal_var[[i]])!=nSNPs[[i]])cat(paste0((nSNPs[[i]]-sum(nSNPs_causal_var[[i]]))," SNPs are assigned zero effects for Trait",i,"! \n"))
        nSNPs_causal[[i]]=nSNPs_causal_var[[i]]
      }
    }
  
  }
    
  return(nSNPs_causal)
  
}

#define the cov matrix of genetic correlation
cal_SNPs_cov<-function(nTraits=2,SNPs_archi=NULL,nSNPs=NULL,h2SNPs=NULL,nSNPs_causal_var=NULL,SNPs_cor=NULL){
  SNPs_cov_mat=NULL
  #for correlated traits, they shared the same SNPs_archi
  if(SNPs_archi[[1]]=="Infinite"|SNPs_archi[[1]]=="Sparse"){
    SNPs_cov_mat=diag(nTraits)
    for(i in 1:(nTraits)){
      SNPs_cov_mat[i,i]=h2SNPs[[i]]/(nSNPs_causal_var[[i]])
      for(j in 1:i){
        if(j!=i)SNPs_cov_mat[i,j]=SNPs_cov_mat[j,i]=SNPs_cor[j]*(sqrt(h2SNPs[[i]]*h2SNPs[[j]]))/nSNPs_causal_var[[i]] # (nSNPs_causal_shared[[1]])
      }
    }  
  }else if(SNPs_archi[[1]]=="BSLMM"|SNPs_archi[[1]]=="BayesR"){
    SNPs_cov_mat=NULL
    for(k in 1:length(SNPs_causal_var[[1]])){ # for large, small groups respectively
      i_SNPs_cov_mat=diag(nTraits)
      
      for(i in 1:nTraits){
        sigma2=SNPs_causal_var[[i]][k] # SNPs_causal_var is a list
        i_SNPs_cov_mat[i,i]=sigma2*h2SNPs[[i]]/nSNPs_causal_var[[1]][k] #number of SNPs with large effects or small effects
        for(j in 1:i){
          #for different groups, we assume 
          if(j!=i)i_SNPs_cov_mat[i,j]=i_SNPs_cov_mat[j,i]=sigma2*SNPs_cor[j]*(sqrt(h2SNPs[[i]]*h2SNPs[[j]]))/nSNPs_causal_var[[1]][k]  # using all causal SNPs in each group instead of (nSNPs_causal_shared[k]) 
        }
      }
      SNPs_cov_mat=c(SNPs_cov_mat,list(i_SNPs_cov_mat))
    }  
  }
   return(SNPs_cov_mat)
}


#sample beta effects 
#currently, we assumed common SNPs are the same for the genetic correlation
beta_sample_SNPs<-function(nTraits=2,SNPs_archi=NULL,nSNPs=NULL,h2SNPs=NULL,nSNPs_causal_var=NULL,
                               SNPs_cor=NULL,nSNPs_causal_shared=NULL,SNPs_list=NULL,seed=1998){
      #sample beta effects 
      beta_traits=matrix("NA",nrow=max(sapply(nSNPs_causal_var,sum)),ncol=nTraits*3) #index，groups，beta
      
      colnames(beta_traits)[3*(1:nTraits)]=paste0("SNPs_effect_Trait",(1:nTraits))
      colnames(beta_traits)[3*(1:nTraits)-1]=paste0("SNPs_arch_",sapply(SNPs_archi,base::c))
      colnames(beta_traits)[3*(1:nTraits)-2]=paste0("SNPs_index_Trait",(1:nTraits))
      
      set.seed(seed)
      if(!is.null(SNPs_cor)){# for genetic correlated traits
      
        SNPs_cov_mat=cal_SNPs_cov(nTraits,SNPs_archi,nSNPs,h2SNPs,nSNPs_causal_var,SNPs_cor)
        
        #SNPs_common_list=Reduce(SNPs_list),we assumed common SNPs are the same for the genetic correlation
        
        #define SNPs_cov_mat for multivariate sampling
        #for correlated traits, they shared the same SNPs_archi
        if(SNPs_archi[[1]]=="Infinite"|SNPs_archi[[1]]=="Sparse"){
  
          #if(SNPs_causal_mode[[i]]=="Random"){ #assume nSNPs_causal_shared is a numeric value, all traits share the same nSNPs_causal_shared
            #sample indexs of SNPs 
            
            
            if(SNPs_causal_mode[[1]]=="Random"){
              beta_traits[1:sum(nSNPs_causal_var[[1]]),3*(1:nTraits)-2]=sample(c(1:nSNPs[[1]]),sum(nSNPs_causal_var[[1]]),replace = F)
            }  
            
            beta_traits[1:nSNPs_causal_shared,3*(1:nTraits)-1]="Shared_group" # groups information
            #sample shared SNPs effects 
            beta_traits[1:nSNPs_causal_shared,3*(1:nTraits)]=MASS::mvrnorm(nSNPs_causal_shared,rep(0,nTraits),Sigma=(SNPs_cov_mat)) 

            #sample non-shared SNPs effects
            for(i in 1:(nTraits)){
              i_nSNPs_causal_var=nSNPs_causal_var[[i]] # number of all causal_SNPs for trait i
              if(nSNPs_causal_shared<i_nSNPs_causal_var){ #not all causal SNPs are shared across traits
                beta_traits[(nSNPs_causal_shared+1):i_nSNPs_causal_var,i*3]=rnorm(i_nSNPs_causal_var-nSNPs_causal_shared,sd=sqrt(h2SNPs[[i]]/i_nSNPs_causal_var))
                beta_traits[(nSNPs_causal_shared+1):i_nSNPs_causal_var,i*3-1]="non_Shared_group"
              }
              }
           # }
          
          
        }else if(SNPs_archi[[1]]=="BSLMM"|SNPs_archi[[1]]=="BayesR"){
  
          #if(SNPs_causal_mode[[1]]=="Random"){ #assume nSNPs_causal_shared is a vector, all traits share the same nSNPs_causal_shared
            k_index=cumsum(nSNPs_causal_shared) #cumulative of nSNPs_causal_shared
            beta_traits[1:nSNPs_causal_shared[1],3*(1:nTraits)]=MASS::mvrnorm(nSNPs_causal_shared[1],rep(0,nTraits),Sigma=(SNPs_cov_mat[[1]]))# k=1
            beta_traits[1:nSNPs_causal_shared[1],3*(1:nTraits)-1]=("Shared_group1") #shared group1
            
            if(SNPs_causal_mode[[1]]=="Random"){
              beta_traits[1:sum(nSNPs_causal_var[[1]]),3*(1:nTraits)-2]=sample(c(1:nSNPs[[1]]),sum(nSNPs_causal_var[[1]]),replace = F)
            } 
            
            for(k in 2:length(nSNPs_causal_shared)){
              sigma2=SNPs_causal_var[[1]][k] # SNPs_causal_var is a list
              #sample shared SNPs effects 
              beta_traits[(k_index[k-1]+1):k_index[k],1:nTraits]=MASS::mvrnorm(nSNPs_causal_shared[k],rep(0,nTraits),Sigma=(SNPs_cov_mat[[k]]))
              beta_traits[(k_index[k-1]+1):k_index[k],nTraits+(1:nTraits)]=paste0("Shared_group",k) #shared group k
            }
              #sample non-shared SNPs effects
              #nSNPs_causal_shared=c(1000,1000,1000);nSNPs_causal_var=c(2000,3000,4000);i_non_shared=nSNPs_causal_var-nSNPs_causal_shared;k_i_non_shared=c(0,cumsum(i_non_shared))+sum(nSNPs_causal_shared)
              for(i in 1:(nTraits)){
                i_non_shared=nSNPs_causal_var[[i]]-nSNPs_causal_shared
                k_i_non_shared=c(0,cumsum(i_non_shared))+sum(nSNPs_causal_shared)
                for(k in 1:length(nSNPs_causal_shared)){
                  i_nSNPs_causal_var=nSNPs_causal_var[[i]][k] # number of all causal_SNPs for trait i on kth group
                  if(nSNPs_causal_shared[k]<i_nSNPs_causal_var){ #not all causal SNPs are shared across traits
                    beta_traits[(k_i_non_shared[k]+1):(k_i_non_shared[k+1]),i*3]=rnorm(i_nSNPs_causal_var-nSNPs_causal_shared[k],sd=sqrt(SNPs_cov_mat[[k]][i,i]))
                    beta_traits[(k_i_non_shared[k]+1):(k_i_non_shared[k+1]),i*3-1]=paste0("non_Shared_group",k) #shared group k
                  }
              }
            }        
         # }
        }
    }else{ #for independent traits
      
      for(i in 1:nTraits){
        
        if(SNPs_archi[[i]]=="Infinite"|SNPs_archi[[i]]=="Sparse"){
          
          #if(SNPs_causal_mode[[i]]=="Random"){ #assume nSNPs_causal_shared is a numeric value, all traits share the same nSNPs_causal_shared
            beta_traits[1:nSNPs_causal_var[[i]],i*3]=rnorm(nSNPs_causal_var[[i]],sd=sqrt(h2SNPs[[i]]/nSNPs_causal_var[[i]]))
            beta_traits[1:nSNPs_causal_var[[i]],i*3-1]="group"
         # }
          
          
        }else if(SNPs_archi[[i]]=="BSLMM"|SNPs_archi[[i]]=="BayesR"){
          sigma2=SNPs_causal_var[[i]][1] # SNPs_causal_var is a list, sigma2 for the 1st group
          #if(SNPs_causal_mode[[i]]=="Random"){ #assume nSNPs_causal_shared is a vector, all traits share the same nSNPs_causal_shared
              k_index=cumsum(nSNPs_causal_var[[i]]) #cumulative of nSNPs_causal_var
              beta_traits[1:k_index[1],i*3]=rnorm(nSNPs_causal_var[[i]][1],sd=sqrt(h2SNPs[[i]]*sigma2/nSNPs_causal_var[[i]][1])) # k=1
              beta_traits[1:k_index[1],i*3-1]=("group1") #shared group1
              for(k in 2:length(nSNPs_causal_var[[i]])){
                sigma2=SNPs_causal_var[[i]][k] # SNPs_causal_var is a list,  sigma2 for the kth group
                #sample  SNPs effects within each group 
                beta_traits[(k_index[k-1]+1):k_index[k],i*3]=rnorm(nSNPs_causal_var[[i]][k],sd=sqrt(h2SNPs[[i]]*sigma2/nSNPs_causal_var[[i]][k]))
                beta_traits[(k_index[k-1]+1):k_index[k],i*3-1]=paste0("group",k) # group k
              }
            #}
          }
        }        
        
    }
      return(beta_traits)
}


#sample environments, errors
sample_residual<-function(nTraits=2,Res_cor=NULL,h2SNPs=NULL,nINDs=NULL,seed=1998){
  
  residual_traits=matrix(NA,nrow=max(sapply(nINDs,sum)),ncol=nTraits)
  set.seed(seed)
  if(!is.null(Res_cor)){ # for same sets of individuals  
    
    Res_cor_mat=diag(nTraits)
    for(i in 1:nTraits){
      Res_cor_mat[i,i]=1-h2SNPs[[i]]
      for(j in 1:i){
        if(j!=i){Res_cor_mat=Res_cor[j]}
      }
    }
    
    residual_traits[1:nINDs[[i]],1:nTraits]=MASS::mvrnorm(nSNPs_causal_shared[1],rep(0,nTraits),Sigma=Res_cor_mat)
    
  }else{ # for different sets of individuals  
    
    for(i in 1:nTraits){
      
      residual_traits[1:nINDs[[i]],i]=rnorm(nINDs[[i]],sd=sqrt(1-h2SNPs[[i]])) #assume no other variance parts
      
    }
    
  }
  return(residual_traits)
}
