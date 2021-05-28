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
library(ggpubr)


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

linearfit <- function(data, ...){
  mod <- lm(y ~., data)
}

linearpred <- function(mod, newdata){
  as.vector(predict(mod, newdata))
}

randomforestfit <- function(data, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = 100, importance = TRUE, ...)
  rf
}

randomforestpredict <- function(data, newdata){
  as.vector(predict(data, newdata=newdata))
}


sim_data <- function(nstudies, ncoef, ntest){
  studies_list  <- vector("list", nstudies)
  nchoose <- 10
  #nchoose <- ncoef
  #general predictor-outcome rule: 
  coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5))) 
  vars <- sample(1:ncoef, nchoose)
  icoefs <- c(4, 1.8)
  n.noise <- 5
  #If linear outcome model:
  #icoefs <- c(0, 0)
  
  #norm.betas <- norm_vec(c(coefs,icoefs))
  #coefs <- coefs/norm.betas
  #cicoefs <- icoefs/norm.betas
  
  m <- genRandomClust(numClust = nstudies - ntest,
                      sepVal=0.8,
                      numNonNoisy=ncoef - n.noise,
                      numNoisy=n.noise,
                      numOutlier=100,
                      numReplicate=1,
                      fileName="test",
                      clustszind=1,
                      clustSizeEq=500,
                      rangeN=c(300,350),
                      clustSizes=NULL,
                      covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
                      rangeVar=c(1, 10))
  
  m2 <- genRandomClust(numClust = 2,
                       sepVal=0.8,
                       numNonNoisy=ncoef - n.noise,
                       numNoisy=n.noise,
                       numOutlier=50,
                       numReplicate=ntest,
                       fileName="test",
                       clustszind=1,
                       clustSizeEq=500,
                       rangeN=c(500,600),
                       clustSizes=NULL,
                       covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
                       rangeVar=c(1, 10))
  
  outliers <- scale(as.data.frame(m$datList)[which(m$memList[[1]] == 0), ])
  outliers_list <- vector("list", nstudies)
  
  d <- 1:nrow(outliers)
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/(nrow(outliers)/nstudies)))
  for (i in 1:length(d1)){
    outliers_list[[i]] <- as.data.frame(outliers[d1[[i]],])
  }
  
  for (i in 1:nstudies){
    curcoefs <- sapply(coefs, function(x){runif(1, x - .5, x + .5)})
    
    if (i <= (nstudies - ntest)){
      data_i <- as.data.frame(m$datList)[which(m$memList[[1]] == i), ]
      data_i  <- scale(rbind(data_i, outliers_list[[i]]))
    }
    else{
      data_i <-scale(as.data.frame(m2$datList[[i - (nstudies - ntest)]]))
    }
    #data_i <- scale(as.data.frame(m$X[m$id == i, ]))
    
    
    #generate outcome
    #scaled here, but the original variables are left unscaled for the clustering step
    y <- as.matrix((data_i[,vars]) %*% curcoefs) + 
      #icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
      cbind(rnorm(nrow(data_i))) # Added noise

    studies_list[[i]] <- as.data.frame(cbind(y, data_i))
    colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
  }
  return(list(studies_list = studies_list))
}

silhouette_score <- function(K_ls, merged){
  ss_ls <- c()
  km_ls <- vector("list", length(K_ls))
  for (i  in 1:length(K_ls)){
    k <- K_ls[i]
    km <- kmeans(merged, centers = k, nstart=25)
    ss <- mean(silhouette(km$cluster, dist(merged))[,3])
    ss_ls <- c(ss_ls, ss)
    km_ls[[i]] <- km$cluster
  }
  index_max <- which.max(ss_ls)
  #print(index_max)
  return(list(km_ind = km_ls[[index_max]], index_max = index_max))
}


create_clusters <- function(studies_list, ntest, ncoef){
  merged <- do.call(rbind, studies_list[1:(length(studies_list) - ntest)])
  merged <- merged[sample(nrow(merged)), ]
  #cluster without using y
  ss <- silhouette_score(2:(floor(nrow(merged)/100)), merged[,-1])
  #k2 <- kmeans(merged[,-1], centers = nclusters, nstart = 25)$cluster
  k2 <- ss$km_ind
  index_max <- ss$index_max
  clusters_list <- lapply(split(seq_along(k2), k2), #split indices by a
                          function(m, ind) m[ind,], m = merged)[order(unique(k2))]
  nclusters <- length(clusters_list)
  for (t in 1:ntest){
    clusters_list[[(nclusters + t)]] <- studies_list[[(length(studies_list) - ntest + t)]]
  }
  return(list(clusters_list = clusters_list, index_max = index_max))
}


create_random <- function(studies_list, ntest, ncoef){
  ntr <- length(studies_list) - ntest
  edat <-  vector("list", length(studies_list))
  merged <- do.call(rbind, studies_list[1:ntr])
  merged <- merged[sample(nrow(merged)), ]
  
  
  d <- (1:nrow(merged))[sample(1:nrow(merged))]
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/(nrow(merged)/ntr)))
  for (i in 1:length(d1)){
    edat[[i]] <- as.data.frame(merged[d1[[i]],])
  }
  for (t in 1:ntest){
    edat[[(ntr + t)]] <- as.data.frame(studies_list[[(length(studies_list) - ntest + t)]])
  }
  return(edat)
}

#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating substudies, 3 if just running on the original set of simulated studies
#rf_ind = 1 if the SSL is random forest, 2 if the SSL is ridge 
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, rf_ind){
  #edat <- sim_data(ndat, ncoef)$studies_list
  if (cluster_ind == 1){ #K-means clustering
    #edat <- create_clusters(studies_list, ndat - ntest, ntest)$clusters_list
    cc <- create_clusters(studies_list, ntest, ncoef)
    edat <- cc$clusters_list
    index_max <- cc$index_max
  }
  else if (cluster_ind == 2){ #random 'pseudo-studies'
    edat <- create_random(studies_list, ntest, ncoef)
  }
  else if (cluster_ind == 3){ #multistudy learning KNOWING the true clusters
    #edat <- studies_list
    edat <- studies_list
  }
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  
  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])

  if (rf_ind == 1){
    mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  }
  else if (rf_ind == 2){
    mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  }
  
  
  ####Inverse variance weights
  invvar <- rep(NA, ntrain)
  for (t in 1:ntrain){
    xt <- as.matrix(edat[[t]][,-1])
    new_obs <- as.matrix(do.call(rbind, edat[as.vector(1:ntrain)[-t]])[,-1])
    preds_var <- sum(diag(new_obs %*% solve(t(xt)%*% xt) %*% t(new_obs)))
    invvar[t] <- preds_var
  }
  coefs_inv <- invvar/sum(invvar)
  
  
  for (j in 1:ntrain){
    mods[[j]] <- modfit(edat[[j]])
    
    preds <- lapply(edat[1:ntrain], function(x){
      modpred(mods[[j]], newdata = x[, -1]) 
    })
    mses[j,] <- unlist(lapply(edat[1:ntrain], function(x){#cross validation within the training set
      newdata = x[, -1]
      preds <- modpred(mods[[j]], newdata = newdata) 
      mean((preds - x[,"y"])^2)}
    ))
    
    curpreds <- lapply(preds, as.numeric)
    allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
  }
  diag(mses) <- NA
  
  # CS Weights
  tt <- apply(mses, 1, mean, na.rm = T) #removing the diagonal elements, takes the mean of each row 
  weights <- absnorm(sqrt(tt), max.norm = TRUE)
  nk <- unlist(lapply(edat, nrow)) #vector of number of rows in each dataset
  nwts <- absnorm(nk[1:ntrain])
  
  # Regression: stacked (intercept and no intercept)
  predstack <- do.call(rbind, allpreds)
  coefs_stack_noint <- nnls::nnls(predstack, as.numeric(as.character(matstack$y)))$x
  coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), as.numeric(as.character(matstack$y)))$x
  coefs_stack_lasso <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 1, lower.limits = 0, intercept = T)))
  coefs_stack_ridge <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 0, lower.limits = 0, intercept = T)))
  #Just a safeguard against full collinearity, although I think we are OK with nnls now
  coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
  coefs_stack_int[which(is.na(coefs_stack_int))] <- 0
  coefs_stack_lasso[which(is.na(coefs_stack_lasso))] <- 0
  coefs_stack_ridge[which(is.na(coefs_stack_ridge))] <- 0
  
  coefs_stack_noint_norm <- absnorm(coefs_stack_noint)
  
  # Regression: study-specific (intercept and no intercept)
  coefs_ss_noint <- mapply(function(x,y){nnls::nnls(y,as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_noint <- colMeans(do.call(rbind, coefs_ss_noint), na.rm = T)
  coefs_ss_int <- mapply(function(x,y){nnls::nnls(cbind(rep(1,nrow(y)),y),as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_int <- colMeans(do.call(rbind, coefs_ss_int), na.rm = T)
  coefs_ss_lasso <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 1, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_ridge <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 0, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_lasso <- colMeans(do.call(rbind, coefs_ss_lasso), na.rm = T)
  coefs_ss_ridge <- colMeans(do.call(rbind, coefs_ss_ridge), na.rm = T)
  coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
  coefs_ss_int[which(is.na(coefs_ss_int))] <- 0
  coefs_ss_lasso[which(is.na(coefs_ss_lasso))] <- 0
  
  coefs_ss_noint_norm <- absnorm(coefs_ss_noint)
  
  
  #########
  #Average difference between inverse variance weights and stacking weights 
  
  diff_inv_stack <- (sapply(1:ntrain, function(i) ((coefs_stack_ridge[-1][i] - coefs_inv[i])/coefs_inv[i])*100))
  diff_inv_unweighted <- (sapply(1:ntrain, function(i) (((1/ntrain) - coefs_inv[i])/coefs_inv[i])*100))
  diff_stack_unweighted <- (sapply(1:ntrain, function(i) (((1/ntrain) - coefs_stack_ridge[-1][i])/coefs_stack_ridge[-1][i])*100))
  
  diff_total <- rbind(diff_inv_stack, diff_inv_unweighted, diff_stack_unweighted)
  rownames(diff_total) <- c("inv_stack", "inv_unweighted", "stack_unweighted")
  
  
  #IF USING MEAN AND SD 
  # diff_total <- data.frame(matrix(NA, 2, 3))
  # colnames(diff_total) <- c("inv_stack", "inv_unweighted", "stack_unweighted")
  # rownames(diff_total) <- c("mean", "sd")
  # 
  # diff_total[1, ] <- c(mean(diff_inv_stack), mean(diff_inv_unweighted), mean(diff_stack_unweighted))
  # diff_total[2, ] <- c(sd(diff_inv_stack), sd(diff_inv_unweighted), sd(diff_stack_unweighted))
  
  outmat <- matrix(NA, ntest, 4)
  colnames(outmat) <- c("Merged", "Unweighted","Inverse_Var","Stack_ridge")
  
  for(i in (ntrain + 1):(length(edat))){
    merged <- modpred(mod0, newdata = edat[[i]][,-1])
    
    merged <- as.vector(sapply(merged, as.numeric))
    allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
    allmod <- apply(allmod, 2, as.numeric)
    # Unweighted average
    unweighted <- colMeans(allmod)
    
    # sample size weighted
    sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})
    
    #inverse variance weighted 
    inv_wtd <- apply(allmod, 2, function(x){sum(coefs_inv*x)})
    
    # cross-study weighted 
    cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})
    
    # regression: stacked (noint, int, each normed) + lasso
    stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
    stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})
    
    # regression: study_specific (noint, int, noint normed) + lasso
    ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
    ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
    ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})
    ss_lasso <- apply(allmod, 2, function(x){coefs_ss_lasso[1] + sum(coefs_ss_lasso[-1]*x)})
    ss_ridge <- apply(allmod, 2, function(x){coefs_ss_ridge[1] + sum(coefs_ss_ridge[-1]*x)})
    
    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    # outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
    #                               mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
    #                               mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), 
    #                               mean((cury - ss_int)^2), mean((cury - stack_lasso)^2), mean((cury - ss_lasso)^2), 
    #                               mean((cury - stack_ridge)^2), mean((cury - ss_ridge)^2)))
    # 
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2),mean((cury - unweighted)^2),
                                  mean((cury - inv_wtd)^2), mean((cury - stack_ridge)^2)))
  }
  #outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  coefs_output <- rbind(coefs_inv, coefs_stack_ridge[-1])
  rownames(coefs_output) <- c("coefs_inv", "coefs_stack_ridge")
  
  if (cluster_ind == 1){
    return(list(outmat = colMeans(outmat), index_max = index_max))
  }
  else{
    return(list(outmat = colMeans(outmat), diff_total = diff_total, norm_stack = sum(coefs_stack_ridge), coefs_output = coefs_output))
  }
  #return(outmat)
}

rep.clusters_fit <- function(reps, modfit, modpred, ndat, ntest, ncoef){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  #num.threads <- as.integer(reps/10)# round up
  num.threads <- 2
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
    sd <- sim_data(nstudies = ndat, ncoef = ncoef, ntest = ntest)$studies_list
    #for linear regression: 
    errors_linear <- clusters_fit(modfit = linearfit, modpred = linearpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, rf_ind = 2)
    
    return(list(errors_linear$outmat, errors_linear$diff_total[1,], errors_linear$diff_total[3,], errors_linear$norm_stack,
                errors_linear$coefs_output[1, ], errors_linear$coefs_output[2, ]))
  }
  closeAllConnections()
  
  errors_linear = na.omit(results[[1]])
  diff_inv_stack = results[[2]]
  diff_stack_unweighted = results[[3]]
  norm_stack = results[[4]]
  coefs_inv = results[[5]]
  coefs_stack = results[[6]]
  
  colnames(errors_linear) <- c("Merged", "Unweighted","Inverse_Var","Stack_ridge")
  
  means_linear <- colMeans(errors_linear)
  
  sds_linear <- apply(errors_linear, 2, sd)
  
  return(list(means_linear = means_linear,
              sds_linear = sds_linear, 
              errors_linear = errors_linear,
              diff_inv_stack = diff_inv_stack, diff_stack_unweighted = diff_stack_unweighted, 
              norm_stack = norm_stack, 
              coefs_inv = coefs_inv, coefs_stack = coefs_stack))   
}

#Run to obtain results
t2 <- rep.clusters_fit(100, modfit = linearfit, modpred = linearpred, ndat = 10, ntest = 5, ncoef = 20)


#To create the table: 
#First, rank the coefficients:
order_inv <- t(sapply(1:100, function(i) t2$coefs_inv[i,][order(t2$coefs_inv[i,])]))
order_stack <- t(sapply(1:100, function(i) t2$coefs_stack[i,][order(t2$coefs_inv[i,])]))
order_diff_inv_stack <- t(sapply(1:100, function(i) t2$diff_inv_stack[i,][order(t2$coefs_inv[i,])]))
order_diff_stack_unweighted <- t(sapply(1:100, function(i) t2$diff_stack_unweighted[i,][order(t2$coefs_inv[i,])]))


median_inv <- round(apply(order_inv, 2, median), 2)
range_inv <- round(apply(order_inv, 2, range), 2)
quantile_inv <- round(sapply(1:5, function(i) quantile(order_inv[,i], c(.25, .75))), 3)

median_stack <- round(apply(order_stack, 2, median), 2)
range_stack <- round(apply(order_stack, 2, range), 2)


median_diff_inv_stack <- round(apply(order_diff_inv_stack, 2, median), 2)
range_diff_inv_stack <- round(apply(order_diff_inv_stack, 2, range), 2)

median_diff_stack_unweighted <- round(apply(order_diff_stack_unweighted, 2, median), 2)
range_diff_stack_unweighted <- round(apply(order_diff_stack_unweighted, 2, range), 2)



summary_coefs_inv <- sapply(1:5, function(i) paste0(median_inv[i], " (", range_inv[1,i], ", ", range_inv[2,i], ")"))
summary_coefs_stack <- sapply(1:5, function(i) paste0(median_stack[i], " (", range_stack[1,i], ", ", range_stack[2,i], ")"))
summary_diff_inv_stack <- sapply(1:5, function(i) paste0(median_diff_inv_stack[i], " (", range_diff_inv_stack[1,i], ", ", range_diff_inv_stack[2,i], ")"))
summary_diff_stack_unweighted <- sapply(1:5, function(i) paste0(median_diff_stack_unweighted[i], " (", range_diff_stack_unweighted[1,i], ", ", range_diff_stack_unweighted[2,i], ")"))


summary_coefs <- rbind(summary_coefs_inv, summary_coefs_stack, summary_diff_inv_stack)
rownames(summary_coefs) <- c("Inverse Variance (IV)", "Stack Ridge (SR)", "% Difference b/w SR and IV")
colnames(summary_coefs) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")


means_linear <- round(t2$means_linear, 3)
sds_linear <- round(t2$sds_linear, 3)

summary_performance <- as.data.frame(t(sapply(1:4, function(i) paste0(means_linear[i], " (", sds_linear[i], ")"))))
colnames(summary_performance) <- c("Merged", "Simple Average", "Inverse Variance", "Stack Ridge")
rownames(summary_performance) <- "Average RMSE"


# Add titles and footnote
# Wrap subtitle into multiple lines using strwrap()
main.title <- "Summary of coefficients"
tab <- ggtexttable(summary_coefs, theme = ttheme("light"))
tab <- tab %>%
  tab_add_title(text = main.title, face = "bold", padding = unit(0.5, "line"))


main.title1 <- "Performance of approaches"
tab1 <- ggtexttable(summary_performance, theme = ttheme("light"))
tab1  <- tab1 %>%
  tab_add_title(text = main.title1, face = "bold", padding = unit(0.5, "line")) 

#Final table
ggarrange(tab1, tab, ncol = 1, nrow = 2, align = 'hv',  labels  = c("A", "B"))

