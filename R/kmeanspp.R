# My own implementation of the WCSS
my_center_within_ss <- function(data, centers) {
  
  # Calculate the distances of each point in data, to each point in centers
  distToCenter <- apply(centers,1,function(center) {rowSums(sweep(data,2,center)^2)})
  
  # Pick only the minimum distance for each point in data, i.e. closest center
  distToClosestCenter <- apply(distToCenter,1,min)
  
  # sum it and return
  sum(distToClosestCenter)
}

# k-means++ function, originally from the LICOR library
# Found a bug in the computation of the "dists", and fixed below
kmeanspp <- function (data, k = 2, start = "random", iter.max = 100, nstart = 10, ...) 
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
  for (restart in seq_len(nstart)) {
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
    tmp.out <- kmeans(data, centers = data[center_ids, ], 
                      iter.max = iter.max, ...)
    tmp.out$inicial.centers <- data[center_ids, ]
    if (tmp.out$tot.withinss < out$tot.withinss) {
      out <- tmp.out
    }
  }
  invisible(out)
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

# This is the original version of the plusplus
# that only concerning the initial centers
# This function is a contrast to the above plusplus function
# for a comparison, to see if the bug is really fixed
ori.plusplus <- function (data, k = 2, start = "random", ...) 
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
                       # rowSums(sweep(data,2,center)^2)
                       rowSums((data - center)^2)
                     })
    }
    else {
      dists <- apply(data[center_ids, ], 1, function(center) {
        # rowSums(sweep(data,2,center)^2)
        rowSums((data - center)^2)
      })
    }
    probs <- apply(dists, 1, min)
    probs[center_ids] <- 0
    center_ids[ii] <- sample.int(num.samples, 1, prob = probs)
  }

  return(data[center_ids, ])
}

# Computing the WCSS of the given cluster, back in its original input space
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

# Functions to generate the random feature by Random Projection
random_features <- function(x, targetD) { 
  randomMatrix <- matrix(rnorm(ncol(x) * targetD, 0, sd=1) * (1/sqrt(targetD)), ncol=targetD)
  out <- x %*% randomMatrix
  
  return(out)
}
