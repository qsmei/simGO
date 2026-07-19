create_cluster_index<-function(p.cluster,n.sample,seed){

  # Number of samples in each cluster
  cluster_sizes <- round(p.cluster * n.sample)

  # Adjust to ensure total matches 1000 due to rounding issues
  cluster_sizes[1] <- cluster_sizes[1] + (n.sample - sum(cluster_sizes))

  # Create a vector of indices for each sample
  all_indices <- seq_len(n.sample)

  # Shuffle the indices randomly
  set.seed(seed) # for reproducibility
  all_indices <- sample(all_indices)

  # Initialize a list to store subset indices
  subset_indices <- list()

  # Generate subset indices for each group
  start <- 1
  for (i in seq_along(cluster_sizes)){
    end <- start + cluster_sizes[i] - 1
    subset_indices[[i]] <- sort(all_indices[start:end])
    start <- end + 1
  }

# subset_indices now contains indices for each group
return(subset_indices)
}




# Function to assign genes within a block with random type
assign_index2blocks <- function(geneIndex, nGenes, min_size, max_size, seed){
  set.seed(seed)
  
  n <- length(geneIndex)
  
  # Check if constraints are feasible
  if (nGenes * min_size > n) {
    stop("Constraints are infeasible: nGenes * min_size exceeds length of geneIndex.")
  }
  if (nGenes * max_size < n) {
    stop("Constraints are infeasible: nGenes * max_size is less than length of geneIndex.")
  }
  
  # Initialize sizes with minimum size for all blocks
  sizes <- rep(min_size, nGenes)
  
  # Calculate the remaining units to distribute
  remaining_units <- n - sum(sizes)
  
  # Calculate the maximum additional units each block can receive
  max_additional <- rep(max_size - min_size, nGenes)
  
  # Distribute the remaining units randomly while respecting the upper bounds
  if (remaining_units > 0) {
    # Indices of blocks that can receive additional units
    blocks_can_receive <- rep(1:nGenes, times = max_additional)
    # Randomly select blocks to receive one unit at a time
    for (i in 1:remaining_units) {
      index <- sample(blocks_can_receive, 1)
      sizes[index] <- sizes[index] + 1
      # Update the blocks_can_receive vector
      if (sizes[index] >= min_size + max_additional[index]) {
        blocks_can_receive <- blocks_can_receive[blocks_can_receive != index]
      }
    }
  }
  
  # Now divide geneIndex into blocks based on sizes
  blocks <- data.frame(matrix(0,nrow=nGenes,ncol=3))
  blocks[,1]=paste0("Gene",1:nGenes)
  start_idx <- 1
  for (i in 1:nGenes) {

    blocks[i,2]=start_idx
    end_idx <- start_idx + sizes[i] - 1
    blocks[i,3]=end_idx
    start_idx <- end_idx + 1

  }
  
  return(blocks)
}



      if(target_dataset=="TCGA_OV"){

        cov.Methy=InterSIM::cov.M
        cov.Gene=InterSIM::cov.expr 
        cov.Protein=InterSIM::cov.protein

        mean.Methy=InterSIM::mean.M
        mean.Gene=InterSIM::mean.expr
        mean.Protein=InterSIM::mean.protein

      }

# generate omics data  based on genotype data (scaled) 
# can also be applied to generate other omics data, if the input data is a omics data
geno2omics<-function(cor_geno_omics=NULL, #the direct correlation between omics1 and omics2
                                           #in theory, the correlation between omic1 and omics2 ~ h2 (e.g.PVEzx)

                     h2_omics=NULL,  #h2 of omics contributed by geno 
                     cor_omics=NULL, # correlation across omics data 
                                     # default is independent, should be identical matrix
                     cov_omics=NULL, # covariance across omics data
                     n.cluster=NULL, #cluster of the number of data
                     p.cluster=c(0.25,0.25,0.25,0.25), #proportions of the sample within each cluster
                     mean.omics.cluster=0,
                     delta.omics.cluster=5,
                     p.deomics=0.2, #proportions of the omics features are differentialy expressed
                                    #if no specific, the p.deomics are the same across p.cluster
                                    #but the index of omics features would be differet across differnt clusters
                     geno=NULL, #omics1 data, should have been scaled 
                     omics_type=NULL, # c("Methy","Gene","Protein")

                     align_block_index=NULL, # dataframe format with 3 columns or 2 columns , the first column is feature in omics2 
                                             # for 3 columns format: 2 and 3 are the start and the end index of the corresponding feature in omics1
                                             # for 2 columns format: 2 is the index of the corresponding feature in omics1, that mean the first column should be repeat

                     align_block_type=NULL,  # random,even,           
                     algin_block_min_size=1, # for random, the minimum size of each bin
                     algin_block_max_size=5, # for random, the minimum size of each bin
                     seed=1998
                     ){

    if(is.null(cov_omics)&is.null(cor_omics))stop("cov_omics or cor_omics shouldn't be NULL!")
    if(is.null(n.cluster)&is.null(p.cluster))stop("p.cluster or n.cluster shouldn't be NULL!")
    if(!is.null(p.cluster))n.cluster=length(p.cluster)
    if(sum(p.cluster)!=1)stop("The sum of p.cluster should equals to 1 !")
    if(is.null(p.cluster)&!is.null(n.cluster))p.cluster=rep(1/n.cluster,n.cluster)
    if(length(p.deomics)==1)p.deomics=rep(p.deomics,n.cluster)

    if(length(delta.omics.cluster)==1)delta.omics.cluster=rep(delta.omics.cluster,n.cluster) #defind the delta shift of mean for each cluster
    if(length(mean.omics.cluster)==1) mean.omics.cluster=rep(mean.omics.cluster,n.cluster)


    if(!is.null(cor_geno_omics)&is.null(h2_omics))h2_omics=cor_geno_omics^2
    #if(!is.null(h2_omics)&is.null(cor_geno_omics))cor_geno_omics=sqrt(cor_geno_omics)

    nIND=nrow(geno) #number of total samples based on the input data-omics1
    nOmics=ncol(cov_omics) #number of Omics features    

    if(is.null(rownames(geno))){
      INDs=paste0("IND",1:nrow(geno))
    }else{
      INDs=rownames(geno)
    }

    group=data.frame(IND=INDs,cluster="C1")

    # index of differenatially expressed CpGs
    set.seed(seed)
    index_deomics <- sapply(1:n.cluster,function(x) rbinom(nOmics, 1, prob = p.deomics[x]))

    cluster_index=create_cluster_index(p.cluster,nIND,seed) #list of sample index corresponding to each cluster

    nSNPs_index=apply(as.matrix(align_block_index[,-1]),1,function(x)x[2]-x[1]+1) #number of snps within each blocks

    # Create an empty list of the specified length
    data_cluster <- vector("list", n.cluster)
    geno_beta_cluster <- vector("list", n.cluster)

    for(i_cluster in 1:n.cluster){

      #the mean effect of omics features across each cluster
      mean.effect <- mean.omics.cluster + index_deomics[,i_cluster]*delta.omics.cluster[i_cluster]

      nIND_cluster=length(cluster_index[[i_cluster]])

      group[cluster_index[[i_cluster]],2]=paste0("C",i_cluster)

      geno_beta=matrix(NA,nrow=nIND_cluster,ncol=nOmics) 
      rownames(geno_beta)=group[cluster_index[[i_cluster]],1]

      for(i in 1:nOmics){

        set.seed(seed+i)
        nSNPs=nSNPs_index[i] #number of SNPs within ith blocks
        i_align_index=align_block_index[i,2]:align_block_index[i,3] #only consider the continuous index 
        subset_geno=as.matrix(geno[cluster_index[[i_cluster]],i_align_index])
        beta <- matrix(rnorm(nSNPs, 0, sd = sqrt(h2_omics[i]/nSNPs)), nSNPs, 1)                                                                                        
        geno_beta[,i]=(subset_geno %*% beta)

      }
      #sample residuals 
      set.seed(seed)
      residual <- mvtnorm::rmvnorm(n = nIND_cluster, mean = rep(0, nOmics), sigma = cov_omics)
      # Combine methyl_data and residual to form the final dataset
      omics <- geno_beta
      for(i in 1:nOmics){
        #omics[,i]=scale(geno_beta[,i]+sqrt(1/h2_omics[i]-1)*sd(geno_beta[,i])*residual[,i]/sd(residual[,i]))+mean.effect[i]
        omics[,i]=scale(geno_beta[,i]+sqrt(1/h2_omics[i]-1)*sd(geno_beta[,i])*residual[,i]/sd(residual[,i]),scale=F)+mean.effect[i] 
      }

      data_cluster[[i_cluster]]=omics

      geno_beta_cluster[[i_cluster]]=geno_beta



      rm(geno_beta,residual,omics);gc(); 

    }

    data_cluster=do.call(rbind,data_cluster)
    geno_beta_cluster=do.call(rbind,geno_beta_cluster)
    # if omics type is "Methylation", rev.logit()
    # if omics type is "Gene expression", generate count data 

    return(list(omics=data_cluster,geno_beta=geno_beta_cluster,group=group,index_deomics=index_deomics))
    }

#example


nIND <- 1000  # Number of individuals
nSNPs <- 367*100  # Number of SNPs, assume 100 SNPs contribute to one CpG 
min_maf <- 0.05  # Minimum MAF threshold
set.seed(123)  # Set seed for reproducibility
mafs <- runif(nSNPs, min = min_maf, max = 0.5)  # MAFs between 0.05 and 0.5



genotype_matrix <- matrix(NA, nrow = nIND, ncol = nSNPs)
# Fill genotype matrix
for (i in 1:nSNPs) {
  p <- mafs[i]  # MAF for current SNP
  set.seed(123+i)
  # Generate genotypes (0, 1, 2) using binomial distribution
  genotype_matrix[, i] <- rbinom(nIND, size = 2, prob = p)
}

setwd("/restricted/projectnb/adiposity/qsmei/Phase3_1000G_imputed/1000GP_Phase3/Hapmap3_1000G/simulation/EUR/Rcpp_code")
Rcpp::sourceCpp("readbin.cpp")
Rcpp::sourceCpp("scaleM.cpp")

genotype_matrix=scaleM_Y(genotype_matrix)


align_block_index=data.frame(ID=paste0("gene",1:367),start=seq(1,10*367,by=10),end=seq(10,10*367,by=10))

result=geno2omics(#cor_geno_omics=NULL, #the direct correlation between omics1 and omics2
                  h2_omics=rep(0.1,367),  #h2 of omics contributed by geno 
                  cov_omics=InterSIM::cov.M, # covariance across omics data
                  mean.omics.cluster=InterSIM::mean.M,
                  n.cluster=4, #cluster of the number of data
                  p.cluster=c(0.25,0.25,0.25,0.25), #proportions of the sample within each cluster
                  delta.omics.cluster=5,
                  align_block_index=align_block_index,
                  p.deomics=0.2, #proportions of the omics features are differentialy expressed
                  geno=genotype_matrix, #omics1 data, should have been scaled 
                  seed=1998
                  )


sim.methyl=result[[1]]

beta_g=result[[2]]
#pos_IND_C1=rownames(sim.methyl)%in%group[group$cluster%in%c("C1"),1]
#sapply(1:ncol(beta_g),function(x)var(beta_g[pos_IND_C1,x])/var(sim.methyl[pos_IND_C1,x]))
group=result[[3]]

index_deomics=result[[4]]

# Load the necessary libraries
library(ggplot2)

# Assume `data` is your original dataset
# Perform PCA
pca_result <- prcomp(sim.methyl, center = TRUE, scale. = TRUE)

# Extract PC1 and PC2 scores for each sample
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

pca_data$cluster=group[match(rownames(sim.methyl),group[,1]),2]
# Plot PC1 vs. PC2 using ggplot
ggplot(pca_data, aes(x = PC1, y = PC2,color=cluster)) +
  geom_point( size = 2, alpha = 0.7) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

library(ggplot2)
library(reshape2)

sim_cov_matrix=cor(sim.methyl[rownames(sim.methyl)%in%group[group$cluster%in%"C2",1],])
#rownames(sim_cov_matrix) <- colnames(sim_cov_matrix) <- paste("CpG", 1:, sep = "_")
rownames(sim_cov_matrix) = colnames(sim_cov_matrix)=NULL
# Convert matrix to long format for ggplot2
cov_melted <- melt(sim_cov_matrix)

# Plot heatmap of correlation
ggplot(cov_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(title = "Simulated: Heatmap of CpG ", x = "CpG", y = "CpG", fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#


library(ComplexHeatmap)


#expr_data_DE <- sim.methyl[, index_deomics[,1]==1 | index_deomics[,2]==1 | index_deomics[,3]==1 | index_deomics[,4]==1]  # Selecting only the first 20 genes
#pos_IND_C1=rownames(sim.methyl)%in%group[group$cluster%in%c("C1"),1]
#postion of non-differentialy expression sites
#pos_non_DE=rowSums(index_deomics) ==0 
#expr_data_non_DE <- sim.methyl[,pos_non_DE] 
#expr_data_DE <- sim.methyl[,!pos_non_DE]  
#rownames(expr_data_non_DE) <- NULL
#colnames(expr_data_non_DE) <- NULL
#rownames(expr_data_DE) <- NULL
#colnames(expr_data_DE) <- NULL
expr_data <- sim.methyl[,] 
rownames(expr_data) <- NULL
colnames(expr_data) <- NULL

Heatmap(expr_data  , name ="Value",column_title = "Heatmap of Simulated Features")
#Heatmap(expr_data_non_DE, name ="Value",column_title = "Heatmap of Non-Differentially Expressed Features")
#Heatmap(expr_data_DE, name ="Value",column_title = "Heatmap of Differentially Expressed Features")




h2_omics=rep(0.1,367)
cov_omics=InterSIM::cov.M # covariance across omics data
n.cluster=4
p.cluster=c(0.25,0.25,0.25,0.25)
mean.omics.cluster=0
delta.omics.cluster=5
align_block_index=align_block_index
p.deomics=0.2
geno=genotype_matrix
seed=1998



h2_SNP2Methy=0.1

pos_sample1=1:200
pos_sample2=201:700 
pos_sample3=701:1000 

pos_sample=list(pos_sample1,pos_sample2,pos_sample3)

delta.methyl=5

methy_effect=mean.M+DMP[,1]*delta.methyl # methy_effect for each CpG within each cluster

#calculat the relationship matrix of methylation
sd_values <- sqrt(diag(cov.M))  # Calculate standard deviations (sqrt of diagonal)
r_methy=cov.M / (sd_values %*% t(sd_values))  # Element-wise division

cov_residual=cov.M


