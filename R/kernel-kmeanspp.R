library(R.matlab)
library(kernlab)
library(tm)

# Function to compute the WCSS back to its input space
sum_dist_to_cluster_center <- function(data, cluster, k) {
  if (nrow(data) != length(cluster)) {
    stop("number of rows in data != length of the given cluster labels");
  }
  
  # Map the result cluster label back to the original space
  # and compute the centroid of the clustered points
  combined <- cbind(data, cluster)
  centersInOrigSpace <- c();
  
  # loop each label
  for (c in 1:k) {
    clusteredPoints <- combined[combined[,ncol(combined)]==c,-ncol(combined)] # find the rows per each label
    centroid <- apply(clusteredPoints, 2, mean)
    centersInOrigSpace <- rbind(centersInOrigSpace, centroid)
  }
  oriSpace_withinss <- my_center_within_ss(data, centersInOrigSpace)
  return(oriSpace_withinss)
}

# This is a light-weighted version of the k-means++
# that only concerning the initial centers. 
# This function would NOT execute the k-means cluster.
plusplus <- function (data, k = 2, start = "random", ...) 
{
  kk <- k
  if (length(dim(data)) == 0) {
    data <- matrix(data, ncol = 1)
  }
  else {
    data <- cbind(data)
  }
  num.samples <- nrow(data)
  ndim <- ncol(data)
  data.avg <- colMeans(data)
  data.cov <- cov(data)
  out <- list()
  out$tot.withinss <- Inf
  
  center_ids <- rep(0, length = kk)
  if (start == "random") {
    center_ids[1:2] = sample.int(num.samples, 1)
  }
  else if (start == "normal") {
    center_ids[1:2] = which.min(dmvnorm(data, mean = data.avg, 
                                        sigma = data.cov))
  }
  else {
    center_ids[1:2] = start
  }
  for (ii in 2:kk) {
    if (ndim == 1) {
      dists <- apply(cbind(data[center_ids, ]), 1, 
                     function(center) {
                       rowSums(sweep(data,2,center)^2)
                       # rowSums((data - center)^2)
                     })
    }
    else {
      dists <- apply(data[center_ids, ], 1, function(center) {
        rowSums(sweep(data,2,center)^2)
        # rowSums((data - center)^2)
      })
    }
    probs <- apply(dists, 1, min)
    probs[center_ids] <- 0
    center_ids[ii] <- sample.int(num.samples, 1, prob = probs)
  }
  
  return(data[center_ids, ])
}

# Objective function of k-means
my_center_within_ss <- function(data, centers) {
  
  # Calculate the distances of each point in data, to each point in centers
  distToCenter <- apply(centers,1,function(center) {rowSums(sweep(data,2,center)^2)})
  
  # Pick only the minimum distance for each point in data, i.e. closest center
  distToClosestCenter <- apply(distToCenter,1,min)
  
  # sum it and return
  sum(distToClosestCenter)
}

random_features <- function(x, targetD) { 
  randomMatrix <- matrix(rnorm(ncol(x) * targetD, 0, sd=1) * (1/sqrt(targetD)), ncol=targetD)
  out <- x %*% randomMatrix
  
  return(out)
}

evalKKmeans.2 <- function(data, k=2, targetDs=c(), sample=0, seed=seed) {
  out <- c()
  
  if (sample>0 && sample < nrow(data)) {
    train.idx <- sample.int(nrow(data), size=sample)
    data <- data[train.idx,]
  }
  
  done <- FALSE
  trial <- 0
  while(!done) {
    trial <- trial+1
    print(paste("Try-Catch loop of trial ", trial, " in samples ", sample, " and seed ", seed, sep=""))
    tryCatch({
      
      # Same seed to make the center picking the same
      ts_pp_rp <- 0
      ts_pp_ori <- 0
      #        set.seed(as.numeric(Sys.time())) # seed by the millis
      
      # Do the plus plus on original dataset with the same seed forumla
      local_seed <- 100*seed + 1*trial
      set.seed(local_seed)
      ts_pp_ori <- system.time(centers.Ori <- plusplus(data=data, k=k))[[1]]
      
      a <- 0
      b <- 0
      x <- 0
      y <- 0
      ts_kkm_rp <- 0
      ts_kkmpp_rp <- 0
      ts_kkm_ori <- 0
      ts_kkmpp_ori <- 0
      ts_randomize <- 0
      
      # Do the kernel kmeans & kernel kemans with plus plus, in original data space
      ts_kkmpp_ori <- system.time(kkmpp.ori <- kkmeans(x=data, centers=centers.Ori, kernel="rbfdot", kpar="automatic", alg="kkmeans"))[[1]]
      b <- sum_dist_to_cluster_center(data, kkmpp.ori@.Data, k)
      
      ts_kkm_ori <- system.time(kkm.ori <- kkmeans(x=data, centers=k, kernel="rbfdot", kpar="automatic", alg="kkmeans"))[[1]]
      a <- sum_dist_to_cluster_center(data, kkm.ori@.Data, k)
      

      out <- rbind(out, t(list("n"=nrow(data), "ori.D"=ncol(data), "targetD"=0, "i"=local_seed, "k"=k,
                               "kkm.withinss"=a, "kkmpp.withinss"=b,
                               "kkm.RP.withinss"=x, "kkmpp.RP.withinss"=y,
                               "kkm.ts"=ts_kkm_ori, "kkmpp.ts"=ts_pp_ori + ts_kkmpp_ori,
                               "kkm.RP.ts"=ts_randomize + ts_kkm_rp, "kkmpp.RP.ts"=ts_randomize + ts_pp_rp + ts_kkmpp_rp, 
                               "kkm.random.ts"=ts_randomize)))
      
      # Loop for the target Ds
      for (D in targetDs) {
        a <- 0
        b <- 0
        x <- 0
        y <- 0
        ts_kkm_rp <- 0
        ts_kkmpp_rp <- 0
        ts_kkm_ori <- 0
        ts_kkmpp_ori <- 0
        ts_randomize <- 0
        ts_pp_rp <- 0
        ts_pp_ori <- 0
        
        print(paste(k,"-",D,sep=""))
        
        # Prepare of random features
        set.seed(as.numeric(Sys.time()))  # need to use a new sequence of random number
        ts_randomize <- system.time(dataRP <- random_features(data, D))[[1]]
        
        # Do the plus plus on projected data 
        set.seed(local_seed)  # need to use the same local seed as the one we used in the ++ on raw input data
        ts_pp_rp <- system.time(centers.RP <- plusplus(data=dataRP, k=k))[[1]]          

        # Do the kernel kmeans & kernel kemans with plus plus, in original data space
        ts_kkmpp_rp <- system.time(kkmpp.rp <- kkmeans(x=dataRP, centers=centers.RP, kernel="rbfdot", kpar="automatic", alg="kkmeans"))[[1]]
        y <- sum_dist_to_cluster_center(data, kkmpp.rp@.Data, k)
        
        ts_kkm_rp <- system.time(kkm.rp <- kkmeans(x=dataRP, centers=k, kernel="rbfdot", kpar="automatic", alg="kkmeans"))[[1]]
        x <- sum_dist_to_cluster_center(data, kkm.rp@.Data, k)
        
        out <- rbind(out, t(list("n"=nrow(data), "ori.D"=ncol(data), "targetD"=D, "i"=local_seed, "k"=k,
                                 "kkm.withinss"=a, "kkmpp.withinss"=b,
                                 "kkm.RP.withinss"=x, "kkmpp.RP.withinss"=y,
                                 "kkm.ts"=ts_kkm_ori, "kkmpp.ts"=ts_pp_ori + ts_kkmpp_ori,
                                 "kkm.RP.ts"=ts_randomize + ts_kkm_rp, "kkmpp.RP.ts"=ts_randomize + ts_pp_rp + ts_kkmpp_rp, 
                                 "kkm.random.ts"=ts_randomize)))
      }
        
      done = TRUE
    }, error = function(err) {
      print(err)
    })
  }
  
  return(out)  
}

out <- c()
data <- readMat("/dataset/GLI_85/GLI_85.mat")$X
loop <- 10
for (k in c(2)) {
  for (i in 1:loop) {
    out <- rbind(out, evalKKmeans.2(data, k=k, targetDs=c(2,3,5,8,10,20,50,100,200,500), sample=0, seed=i))
    write.csv(out, "gli.full.csv")
  }
}

out <- c()
data <- readMat("/dataset/SMK_187/SMK_CAN_187.mat")$X
loop <- 10
for (k in c(2)) {
  for (i in 1:loop) {
    out <- rbind(out, evalKKmeans.2(data, k=k, targetDs=c(2,3,5,8,10,20,50,100,200,500), sample=0, seed=i))
    write.csv(out, "smk.full.csv")
  }
}