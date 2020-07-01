library(kernlab)

## S4 method for signature 'matrix'
# Kernel k-means++ function
kkmeanspp <- function(x, centers, nstart=1, kpar=NULL, ...) {
  kk <- centers
  data <- x

  num.samples <- nrow(data)
  ndim <- ncol(data)

  if (ndim < 2) {
    stop("'x' must be an array of at least two dimensions") 
  }
  
  out <- NULL
  for (restart in seq_len(nstart)) {
    center_ids <- rep(0, length = kk)
    
    # Initialize the first two center
    center_ids[1:2] = sample.int(num.samples, 1)

    # for each center from 2 to k
    for (ii in 2:kk) {
      dists <- apply(data[center_ids, ], 1, function(center) {
         rowSums(sweep(data,2,center)^2)
        # sweep should be equivalent to this
        # rowSums(t(t(data)-center)^2)
      })

      probs <- apply(dists, 1, min)
      probs[center_ids] <- 0   # set probability to 0 for those picked centers
      
      center_ids[ii] <- sample.int(num.samples, 1, prob = probs)
      
      # Detect if any 2 centers are duplicated, redo the sampling using the same probability set
      # This line is to preventing system breakdown when duplicated initial points were chosen
      while (any(duplicated(data[center_ids, ]))) {
        center_ids[ii] <- sample.int(num.samples, 1, prob = probs)
      }
      
    }  # end for loop of centers ii
    
    tmp.out <- kkmeans(data, centers = data[center_ids, ], kpar=kpar, ...)
    
    if (is.null(out) || sum(tmp.out@withinss) < sum(out@withinss)) {
      out <- tmp.out
    }
  } # end for loop of nstart
    
  out
}

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
