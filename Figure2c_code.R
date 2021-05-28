library(locfit)
library(fungible)
library(mvtnorm)
library(ClusterR)
library(cluster)
library(MLmetrics)
library(randomForest)
library(glmnet)
library(nnet)
library(nnls)
library(doParallel)
library(foreach)
library(plyr)
library(BiocGenerics)
library(clusterGeneration)
library(LaplacesDemon)


rcomb <- function(...) {
  args <- list(...)
  lapply(seq_along(args[[1]]), function(i)
    do.call('rbind', lapply(args, function(a) a[[i]])))
}


norm_vec <- function(x) sqrt(sum(x^2))

absnorm <- function(vec, max.norm = FALSE){
  sgn <- sign(vec)
  vec <- abs(vec)
  if(max.norm){
    mvec <- max(vec)	
  } else {
    mvec <- min(vec)
  }
  
  vec <- abs(vec - mvec)
  sgn*(vec/sum(vec))
}

randomforestfit <- function(data, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = 100, importance = TRUE, ...)
  rf
}

randomforestpredict <- function(data, newdata){
  as.vector(predict(data, newdata=newdata))
}



######## Simulating Laplacian data for the Merged
sim_data_merged <- function(nstudies, ncoef, ntest, sampsize){
  studies_list  <-  vector("list", nstudies)
  nchoose <- ncoef
  #nchoose <- 10
  #coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5))) 
  #coefs <- coefs/norm_vec(coefs)
  ntrain <- nstudies - ntest
  
  cluster_means <- t(sapply(1:ntrain, function(i) {
    sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
  }
  ))
  
  coefs <- rnorm(nchoose)
  
  vars <- sample(1:ncoef, nchoose)
  icoefs <- c(4, 1.8)
  numsample <- round(3000/ntrain)
  
  for (i in 1:nstudies){
    curcoefs <- sapply(coefs, function(x){runif(1, x - .25, x + .25)})
    if  (i == nstudies){ #test set
      data_i <- as.matrix(do.call(rbind, lapply(1:ntrain, 
                                                function(t) {rmvl(n = numsample, mu = cluster_means[t, ], Sigma = diag(nchoose))
                                                })))
      
    }
    else{ #training set
      data_i <- rmvl(n = numsample, mu = cluster_means[i, ], Sigma = diag(nchoose))
    }
    
    y <- as.matrix(as.matrix(as.matrix(data_i)[,vars]) %*% curcoefs) + 
      #icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
      cbind(rnorm(nrow(as.matrix(data_i))))
    studies_list[[i]] <- as.data.frame(cbind(y, data_i))
    
    
    colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
  }
  return(list(studies_list = studies_list))
}


######## Simulating Laplacian data for the clusters
sim_data_cluster <- function(nstudies, ncoef, ntest, sampsize){
  studies_list  <-  vector("list", nstudies)
  nchoose <- ncoef

  coefs <- rnorm(nchoose)
  vars <- sample(1:ncoef, nchoose)
  icoefs <- c(4, 1.8)
  
  ntrain <- nstudies - ntest
  
  cluster_means <- t(sapply(1:nstudies, function(i) {
    sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
  }
  ))

  numsample <- round(3000/ntrain)
  
  for (i in 1:nstudies){
    curcoefs <- sapply(coefs, function(x){runif(1, x - .25, x + .25)})
    #curcoefs <- coefs
    if (i == nstudies){
      data_i <- mvrnorm(n = numsample, mu = cluster_means[i, ], Sigma = diag(nchoose))
    }
    else{
      data_i <- rmvl(n = numsample, mu = cluster_means[1, ], Sigma = diag(nchoose))
    }
    y <- as.matrix(as.matrix(as.matrix(data_i)[,vars]) %*% curcoefs) + 
      #icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
      cbind(rnorm(nrow(as.matrix(data_i))))
    studies_list[[i]] <- as.data.frame(cbind(y, data_i))
    
    
    
    colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
  }
  return(list(studies_list = studies_list))
}

#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating substudies, 3 if just running on the original set of simulated studies
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, nsep){
  if (cluster_ind == 3){ #multistudy learning KNOWING the true clusters
    edat <- studies_list
  }
  
  ntrain = length(edat) - ntest
  
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  
  #For centered forest
  mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  
  #For normal forest 
  #mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  outmat <- matrix(NA, ntest, 1)
  colnames(outmat) <- c("Merged")
  
  for(i in (ntrain + 1):(length(edat))){
    #merged <- modpred(mod0, newdata = edat[[i]][,-1])
    merged <- modpred(mod0, newdata = edat[[i]])
    merged <- as.vector(sapply(merged, as.numeric))
    
    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2)))
  }
  return(list(outmat = colMeans(outmat)))
}


rep.clusters_fit <- function(reps, modfit, modpred, ndat, ntest, ncoef, sampsize){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  #num.threads <- as.integer(reps/reps)# round up
  #num.threads <- as.integer(10)
  num.threads <- 15
  threads <- makeCluster(num.threads, outfile=logfile, setup_timeout = 0.5)
  registerDoParallel(threads)
  
  getDoParWorkers()
  #######
  
  results <- foreach(i = 1:reps, .combine = 'rcomb', .multicombine = TRUE, .export = ls(globalenv())) %dopar% {
    library(locfit)
    library(fungible)
    library(mvtnorm)
    library(ClusterR)
    library(cluster)
    library(MLmetrics)
    library(randomForest)
    library(glmnet)
    library(nnet)
    library(nnls)
    library(plyr)
    library(BiocGenerics)
    library(clusterGeneration)
    library(LaplacesDemon)
    sdc <- sim_data_cluster(nstudies = ndat, ncoef = ncoef, ntest = ntest, sampsize = sampsize)$studies_list
    sdm <- sim_data_merged(nstudies = ndat, ncoef = ncoef, ntest = ntest, sampsize = sampsize)$studies_list
    xnam <- paste0("V", 1:ncoef)
    formula <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
    errors_multi_cluster <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sdc, cluster_ind = 3, nsep = 0)$outmat
    errors_multi_merged <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sdm, cluster_ind = 3, nsep = 0)$outmat
    print(errors_multi_cluster)
    print(errors_multi_merged)
    return(list(errors_multi_cluster = errors_multi_cluster, errors_multi_merged = errors_multi_merged))
  }
  closeAllConnections()
  
  errors_multi_cluster = results[[1]]
  errors_multi_merged = results[[2]]
  
  colnames(errors_multi_cluster) <- colnames(errors_multi_merged) <- c("Merged")
  
  means_cluster <- colMeans(errors_multi_cluster)
  sds_cluster <- apply(errors_multi_cluster, 2, sd)
  
  means_merged <- colMeans(errors_multi_merged)
  sds_merged <- apply(errors_multi_merged, 2, sd)
  
  return(list(means_cluster = means_cluster, sds_cluster = sds_cluster, errors_multi_cluster = errors_multi_cluster,
              means_merged = means_merged, sds_merged = sds_merged, errors_multi_merged = errors_multi_merged))   
}


vary_levels <- function(reps, var_list, modfit, modpred, ncoef, ntest, out_str){
  ptm = proc.time()
  colnames_total <- c("Merged")
  total_means_merged <- total_sds_merged <- array(0, c(length(var_list), 1))
  colnames(total_means_merged) <- colnames(total_sds_merged) <-  colnames_total 
  
  total_means_cluster <- total_sds_cluster <- array(0, c(length(var_list), 1))
  colnames(total_means_cluster) <- colnames(total_sds_cluster) <-  colnames_total 
  
  
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    
    level_rep <- rep.clusters_fit(reps, modfit, modpred, ndat = level, ntest, ncoef, sampsize = 500)
    print(level_rep)
    #means
    total_means_merged[i,] <- level_rep$means_merged
    total_means_cluster[i,] <- level_rep$means_cluster
    #sds
    total_sds_merged[i,] <- level_rep$sds_merged
    total_sds_cluster[i,] <- level_rep$sds_merged
    #indices
    write.table(level_rep$errors_multi_merged, paste0(out_str,"_errors_merged",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_multi_cluster, paste0(out_str,"_errors_cluster",level,".csv"), sep = ",", col.names = colnames_total)
  }
  total_means <- cbind(total_means_merged, total_means_cluster)
  total_sds <- cbind(total_sds_merged, total_sds_cluster)
  colnames(total_means) <- colnames(total_sds) <- c("Merged", 'Ensemble')
  write.table(total_means, paste0(out_str,"_means1.csv"), sep = ",", row.names = var_list, col.names = c("Merged", 'Ensemble'))
  write.table(total_sds, paste0(out_str,"_sds1.csv"), sep = ",", row.names = var_list, col.names = c("Merged", 'Ensemble'))
  
  return(list(total_means = total_means, total_sds_multi = total_sds))
}

#Run this to obtain the results
nclust_list <- c(3, 5, 9, 17, 33)
cs1 <- vary_levels(reps = 100, var_list = nclust_list, modfit = randomforestfit, modpred = randomforestpredict, ncoef = 20, ntest = 1, out_str = "cs")


