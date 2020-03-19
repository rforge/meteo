near.obs <- function(
  locations,
  observations,
  zcol = 1,
  n.obs = 50,
  avg = FALSE,
  direct = FALSE,
  increment,
  range,
  rm.dupl = TRUE
)
{
  if (class(locations) == "SpatialPoints" ||
      class(locations) == "SpatialPointsDataFrame" ||
      class(locations) == "SpatialPixelsDataFrame") {
    locations <- coordinates(locations)
  }
  if (class(observations) == "SpatialPoints" || class(observations) == "SpatialPointsDataFrame") {
    variable <- observations[[zcol]]
    observations <- coordinates(observations)
  } else {
    variable <- observations[, zcol]
    observations <- observations[, 1:2]
  }
  
  # if (identical(locations, observations)){
  if (rm.dupl){
    knn1 <- nabor::knn(observations, locations, k=n.obs+1)
    knn1$nn.idx[round(knn1$nn.dists[, 1]) == 0, 1:n.obs] <- knn1$nn.idx[round(knn1$nn.dists[, 1]) == 0, -1]
    knn1$nn.idx <- knn1$nn.idx[, -(n.obs+1)]
    knn1$nn.dists[round(knn1$nn.dists[, 1]) == 0, 1:n.obs] <- knn1$nn.dists[round(knn1$nn.dists[, 1]) == 0, -1]
    knn1$nn.dists <- knn1$nn.dists[, -(n.obs+1)]
  } else {
    knn1 <- nabor::knn(observations, locations, k=n.obs)
  }
  
  near_o1 <- apply(knn1$nn.idx, 2, function(x) {variable[x]})
  near_o1 <- cbind(near_o1)
  
  nl_df <- cbind(knn1$nn.dists, near_o1)
  
  name1 <- c()
  name2 <- c()
  for (i in 1:n.obs) {
    name1 <- c(name1, paste("dist", i, sep = ""))
    name2 <- c(name2, paste("obs", i, sep = ""))
  }
  all_names <- c(name1, name2)
  
  ### AVG - BAND ###
  
  if(avg) {
    
    knn2 <- nabor::knn(observations, locations, k = nrow(observations), radius = range)
    if (rm.dupl){
      knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.idx <- knn2$nn.idx[, -(nrow(observations))]
      knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.dists <- knn2$nn.dists[, -(nrow(observations))]
    }
    
    # knn2$nn.idx <- knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, -1]
    knn2$nn.idx <- ifelse(knn2$nn.idx==0,NA,knn2$nn.idx)
    
    near_o2 <- apply(knn2$nn.idx, 2, function(x) {variable[x]})
    near_o2 <- cbind(near_o2)
    
    nc <- ceiling(range/increment)
    avg_st <- matrix(nrow = nrow(locations), ncol = nc)
    
    for (inc in 1:nc) {
      avg_st[, inc] <- apply(ifelse((knn2$nn.dists <= inc * increment) & (knn2$nn.dists >= (inc-1) * increment), near_o2, NA), 1, "mean", na.rm = T)
      # avg_st[, inc] <- apply(ifelse((knn2$nn.dists <= inc * increment), near_o2, NA), 1, "mean", na.rm = T)
    }
    
    avg_st[, nc] <- ifelse(is.na(avg_st[, nc]), mean(variable), avg_st[, nc])
    for (inc in (nc-1):1){
      avg_st[, inc] <- ifelse(is.na(avg_st[, inc]), avg_st[, (inc+1)], avg_st[, inc])
    }
    
    nl_df <- cbind(nl_df, avg_st)
    name3 <- paste("avg", seq(increment, range + increment-1, increment), sep="")
    all_names <- c(all_names, name3)
    
  }
  
  if(direct) {
    
    knn2 <- nabor::knn(observations, locations, k = nrow(observations))
    if (rm.dupl){
      knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.idx <- knn2$nn.idx[, -(nrow(observations))]
      knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.dists <- knn2$nn.dists[, -(nrow(observations))]
    }
    
    # knn2$nn.idx <- knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, -1]
    # knn2$nn.idx <- ifelse(knn2$nn.idx==0,NA,knn2$nn.idx)
    
    near_o2 <- apply(knn2$nn.idx, 2, function(x) {variable[x]})
    near_o2 <- cbind(near_o2)
    
    # direct_m <- matrix(nrow = nrow(locations), ncol = dim(near_o2)[2])
    direct_obs <- matrix(nrow = nrow(locations), ncol = 4)
    
    direct_obs[, 1] <- sapply(1:nrow(locations),
                                  function(x) near_o2[x, (observations[knn2$nn.idx[x, ], "x"]>locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]>locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]

    direct_obs[, 2] <- sapply(1:nrow(locations),
                                   function(x) near_o2[x, (observations[knn2$nn.idx[x, ], "x"]>locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]<locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]

    direct_obs[, 3] <- sapply(1:nrow(locations),
                                   function(x) near_o2[x, (observations[knn2$nn.idx[x, ], "x"]<locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]<locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]

    direct_obs[, 4] <- sapply(1:nrow(locations),
                                   function(x) near_o2[x, (observations[knn2$nn.idx[x, ], "x"]<locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]>locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    direct_obs <- ifelse(is.na(direct_obs), 999999, direct_obs)
      # avg_st[, inc] <- apply(ifelse((knn2$nn.dists <= inc * increment), near_o2, NA), 1, "mean", na.rm = T)

    nl_df <- cbind(nl_df, direct_obs)
    name3 <- paste("dir", seq(1, 4, 1), sep="")
    all_names <- c(all_names, name3)
    
    
    direct_dist <- matrix(nrow = nrow(locations), ncol = 4)
    
    direct_dist[, 1] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], "x"]>locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]>locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    
    direct_dist[, 2] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], "x"]>locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]<locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    
    direct_dist[, 3] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], "x"]<locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]<locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    
    direct_dist[, 4] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], "x"]<locations[x, "x"] & observations[knn2$nn.idx[x, ], "y"]>locations[x, "y"])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    direct_dist <- ifelse(is.na(direct_dist), 999999, direct_dist)
    
    nl_df <- cbind(nl_df, direct_dist)
    name3 <- paste("dir_dist", seq(1, 4, 1), sep="")
    all_names <- c(all_names, name3)
    
  }
  
  nl_df <- as.data.frame(nl_df)
  names(nl_df) <- all_names
  
  return(nl_df)
  
}









