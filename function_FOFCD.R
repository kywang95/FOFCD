###############################################################################

## All of the procedure code necessary for the paper

###############################################################################

# FunNCC -------------------------------------------------------------------------
#' Estimate the functional conditional dependence coefficient (FunNCC)
FunNCC <- function(Y, Z, X = NULL, na.rm = TRUE, 
                  method=c('fourier_basis', 'ospline_basis', 'raw', 'FPCA'),
                  control=list()){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    if(!is.list(Z)) stop("Z should be the list with p variables, and each variable is a curve.")
    if (sum(diff(sapply(Z, nrow)))!=0) stop("Number of rows of Z should be equal.")
    if (length(Y) != nrow(Z[[1]])) stop("Number of rows of Y and X should be equal.")
    if (na.rm == TRUE) {
      # NAs are removed here:
      ok = complete.cases(Y, Z)
      Z <- lapply(Z, function(x){x[ok,]})
      Y = Y[ok]
    }
    
    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    q = list(Z)
    method = method[1]
    
    obs_Z <- basis_expansion(Z, method, control)
    return(.estimateT(Y, matrix(unlist(obs_Z$obs_X), nrow = n)))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.list(X)) stop("X should be the list with p variables, and each variable is a curve.")
  if(!is.list(Z)) stop("Z should be the list with p variables, and each variable is a curve.")
  rowZ <- rep(0, times = length(Z))
  for(i in 1:length(Z)) rowZ[i] <- nrow(Z[[i]]) 
  if (sum(diff(rowZ))!=0) stop("Number of rows of Z should be equal.")
  rowX <- rep(0, times = length(X))
  for(i in 1:length(X)) rowX[i] <- nrow(X[[i]])
  if (sum(diff(rowX))!=0) stop("Number of rows of X should be equal.")
  if((length(Y) != nrow(X[[1]]))) stop("Number of rows of Y and X should be equal.")
  if((length(Y) != nrow(Z[[1]]))) stop("Number of rows of Y and Z should be equal.")
  if((nrow(Z[[1]]) != nrow(X[[1]]))) stop("Number of rows of Z and X should be equal.")
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(Y,Z,X)
    Z <- lapply(Z, function(x){x[ok,]})
    X <- lapply(X, function(x){x[ok,]})
    Y = Y[ok]
  }
  
  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  q = length(Z)
  p = length(X)
  
  W <- c(X, Z)
  names(W) <- c(1:length(W))
  obs_W <- basis_expansion(W, method, control)
  obs_X <- matrix(unlist(obs_W$obs_X[1:length(X)]), nrow = n)
  obs_Z <- matrix(unlist(obs_W$obs_X[(length(X)+1):length(W)]), nrow = n)
  return(.estimateConditionalT(Y, obs_Z, obs_X))
}



#############################################
# HELPER FUNCTIONS for FunNCC
#
# .estimateConditionalQ
# .estimateConditionalS
# .estimateConditionalT
# .estimateQ
# .estimateS
# .estimateT
##############################################


.estimateConditionalQ <- function (Y, obs_X, obs_Z) {
  
  id <- group <- rnn <- NULL
  
  n = length(Y)
  
  obs_W = cbind(obs_X, obs_Z)
  
  # compute the nearest neighbor of X
  ######################################### RANN::nn2
  nn_X = RANN::nn2(obs_X, query = obs_X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  # handling repeated data
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  
  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      ############################################## Todo: RANN::nn2
      knna <- RANN::nn2(matrix(obs_X[-a,], ncol=ncol(obs_X)), 
                        query = matrix(obs_X[a,], ncol=ncol(obs_X)), 
                        k = nrow(obs_X)-1)
      distances <- knna$nn.idx
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  # compute the nearest neighbor of W
  ######################################### Todo: RANN::nn2
  nn_W = RANN::nn2(obs_W, query = obs_W, k = 3)
  nn_index_W = nn_W$nn.idx[, 2]
  repeat_data = which(nn_W$nn.dists[, 2] == 0)
  
  df_W = data.table::data.table(id = repeat_data, group = nn_W$nn.idx[repeat_data])
  df_W[, rnn := .randomNN(id), by = "group"]
  
  nn_index_W[repeat_data] = df_W$rnn
  # nearest neighbors with ties
  ties = which(nn_W$nn.dists[, 2] == nn_W$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      ############################################## Todo: RANN::nn2
      knna <- RANN::nn2(matrix(obs_W[-a,], ncol=ncol(obs_W)), 
                        query = matrix(obs_W[a,], ncol=ncol(obs_W)), 
                        k = nrow(obs_W)-1)
      distances <- knna$nn.idx
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      # since the a-th observation is removed from W[-a,], if the index x>=a, it should add 1;
      # else, it is the same as the original index.
      return(x + (x >= a))
    }
    
    nn_index_W[ties] = sapply(ties, helper_ties)
  }
  
  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_W]), 2, min),
            -apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  return(Q_n)
}





# .estimateConditionalS -------------------------------------------------------------------------
# Estimate S(Y, X)
#
# Estimate S(Y, X), the denuminator of the measure of dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
#
# @return estimation \eqn{S(Y, X)}
.estimateConditionalS <- function (Y, obs_X){
  
  id <- group <- rnn <- NULL
  
  n = length(Y)
  
  # compute the nearest neighbor of X
  ######################################### Todo: RANN::nn2
  nn_X = RANN::nn2(obs_X, query = obs_X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  
  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      ############################################## Todo: RANN::nn2
      knna <- RANN::nn2(matrix(obs_X[-a,], ncol=ncol(obs_X)),
                        query = matrix(obs_X[a,], ncol=ncol(obs_X)),
                        k = nrow(obs_X)-1)
      distances <- knna$nn.dists
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  # estimate S
  R_Y = rank(Y, ties.method = "max")
  S_n = sum(R_Y - apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  
  return(S_n)
}


# estimateConditionalT -------------------------------------------------------------------------
# Estimate T(Y, Z | X)
#
# Estimate T(Y, Z | X), the measure of dependence of Y on Z given X
#
# @param Y: Vector (length n)
# @param Z: Matrix of predictors (n by q)
# @param X: Matrix of predictors (n by p)
#
# @return estimation of \eqn{T(Y, Z|X)}.
.estimateConditionalT <- function(Y, obs_Z, obs_X){
  S = .estimateConditionalS(Y, obs_X)
  
  # happens only if Y is constant
  if (S == 0) {
    return(1)
  } else {
    return(.estimateConditionalQ(Y, obs_X, obs_Z) / S)
  }
}




# .estimateQ -------------------------------------------------------------------------
# Estimate Q(Y, X)
#
# Estimate Q(Y, X), the numinator of the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p).
# @param Y: Vector (length n).
#
# @return estimation of \eqn{Q(Y, X)}.
.estimateQ <- function(Y, obs_X) {
  
  id <- group <- rnn <- NULL
  
  n = length(Y)
  
  
  ######################################### Todo: RANN::nn2
  nn_X = RANN::nn2(obs_X, query = obs_X, k = 3)
  # remove the first nearest neighbor for each x which is x itself in case of no repeat data
  # when there is repeated data this is wrong but for repeated data we find the nearest
  # neighbors separately.
  nn_index_X = nn_X$nn.idx[, 2]
  
  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  # for the repeated data points, choose one of their identicals at random and set its index
  # as the index of the nearest neighbor
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  nn_index_X[repeat_data] = df_X$rnn
  
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      ############################################## Todo: RANN::nn2
      #      print(obs_X)
      distances <- proxy::dist(matrix(obs_X[-a,], ncol=ncol(obs_X)),
                               matrix(obs_X[a,], ncol=ncol(obs_X)))
      # knna <- RANN::nn2(matrix(obs_X[-a,], ncol=ncol(obs_X)),
      #                   query = matrix(obs_X[a,], ncol=ncol(obs_X)),
      #                   k = nrow(obs_X)-1)
      # distances <- knna$nn.dists
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  
  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  L_Y = rank(-Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_X]), 2, min) - (L_Y^2)/n) / (n^2)
  
  return(Q_n)
}



# .estimateS -------------------------------------------------------------------------
# Estimate S(Y)
#
# Estimate S(Y) , the denuminator of the measure of dependence of Y on X
#
# @param Y: Vector (length n).
# @return estimation of \eqn{S(Y)}.
.estimateS <- function (Y) {
  n = length(Y)
  L_Y = rank(-Y, ties.method = "max")
  S_n = gmp::asNumeric(sum(gmp::as.bigz(L_Y) * gmp::as.bigz(n - L_Y))) / (n^3)
  return(S_n)
}



# .estimateT -------------------------------------------------------------------------
# Estimate T(Y, X)
#
# Estimate T(Y, X), the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
# @return estimation of \eqn{T(Y, X) = Q(Y, X) / S(Y)}.
.estimateT <- function(Y, obs_X){
  # May 15, Mona removed FOCI:::
  S = .estimateS(Y)
  # happens only if Y is a constant vector.
  if (S == 0) {
    return(1)
  } else {
    return(.estimateQ(Y, obs_X) / S)
  }
}


#####################################
# HELPER FUNCTIONS:
# randomNNs
#####################################

# randomNN -------------------------------------------------------------------------
# Find the random nearest neighbors.
#
# Gives the indices of nearest neighbors of points that are not unique in the data set.
# For each point we know that to which cluster of points it belongs (repeated points),
# we chose one of those indices that is not equal to the same index of our given point
# at random as its nearest neighbor.
#
# @param ids: a vector of ids of groups that each index is a member of.
#
# @return a vector of indices.
.randomNN <- function(ids) {
  m <- length(ids)
  
  x <- sample(x = (m - 1), m, replace = TRUE)
  # To avoid ids itself being selected
  x <- x + (x >= (1:m))
  
  return(ids[x])
}


####################################################################
# MAIN FUNCTIONS:
# fofcd_main: the core function for implementing the fofcd algorithm
# fofcd:  Performs feature ordering by functional conditional independence
####################################################################

# fofcd_main -------------------------------------------------------------------------
# Variable selection by the FOCI algorithm
#
# fofcd_main is the core function for foci, implementing the
# variable selection algorithm based on the measure of conditional dependence \code{\link{FunNCC}}.
#
fofcd_main <- function(Y, obs_X, num_features = NULL, stop = TRUE, 
                       numCores = parallel::detectCores(), 
                       method=c('fourier_basis','FPCA','raw'), control=list()){
  
  method = method[1]
  namesX <- names(obs_X)
  if (is.null(num_features)) num_features = length(obs_X)
  n = length(Y)
  p = length(obs_X)
  ndimX = sapply(obs_X, ncol)
  idxFD = seq_along(obs_X)[ndimX > 1]
  
  # obs_X <- basis_expansion(X, method, control)
  
  Q = rep(0, num_features)
  index_select = rep(0, num_features)
  
  # select the first variable
  if (p==1) {
    seq_Q = .estimateQ(Y, obs_X[[1]])
  } else {
    estimateQFixedY <- function(id){
      return(.estimateQ(Y, obs_X[[id]]))
    }
    seq_Q = parallel::mclapply(seq(1, p), estimateQFixedY, mc.cores = numCores)
    seq_Q = unlist(seq_Q)
  }
  
  Q[1] = max(seq_Q)
  if (Q[1] <= 0 & stop == TRUE) return(NULL)
  index_max = min(which(seq_Q == Q[1]))
  index_select[1] = index_max
  
  print(namesX[index_select[1]])
  print(Q[1])
  
  count = 1
  
  # select rest of the variables
  while (count < num_features) {
    seq_Q = rep(0, p - count)
    # indices that have not been selected yet
    index_left = setdiff(seq(1, p), index_select[1:count])
    
    # find the next best feature
    estimateQFixedYandSubX <- function(id){
      obs_W <- matrix(unlist(obs_X[c(index_select[1:count], id)]), nrow = n)
      return(.estimateQ(Y, obs_W))
    }
    
    if (length(index_left) == 1) {
      seq_Q = estimateQFixedYandSubX(index_left[1])
    } else {
      seq_Q = parallel::mclapply(index_left, estimateQFixedYandSubX, mc.cores = numCores)
      seq_Q = unlist(seq_Q)
    }
    Q[count + 1] = max(seq_Q)
    index_max = min(which(seq_Q == Q[count + 1]))
    
    if (Q[count + 1] <= Q[count] & stop == TRUE) break
    index_select[count + 1] = index_left[index_max]
    
    print(namesX[index_select[count + 1]])
    print(Q[count + 1])
    
    count = count + 1
  }
  
  selectedVar = data.table::data.table(index = index_select[1:count], names = namesX[index_select[1:count]])
  stepT = Q / .estimateS(Y)
  result = list(selectedVar = selectedVar, stepT = stepT[1:count])
  class(result) = "foci"
  return(result)
}


# fofcd -------------------------------------------------------------------------
#' Variable selection by the fofcd algorithm
#'
#' FOCI is a variable selection algorithm based on the measure of conditional dependence \code{\link{FunNCC}}.
fofcd <- function(Y, X, num_features = NULL, stop = TRUE, na.rm = TRUE,
                  standardize = "scale", numCores = parallel::detectCores(),
                  parPlat = 'none', printIntermed = TRUE, 
                  method=c('fourier_basis','FPCA','raw'), control=list()) {
  core_time <- Sys.time()
  if(is.null(X)) {
    return(list(feature = NA, time = NA))
  }else{
    if (is.null(names(X))) {
      names(X) <- paste0('x',1:length(X))
      warning('X lacked column names, has been assigned x1, x2,...')
    }
    namesX <- names(X)
    
    # if inputs are not in proper format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y <- as.vector(unlist(Y))
    }
    if(!is.list(X)) stop("X should be the list with p variables, and each variable is a curve.")
    
    if (is.null(num_features)) num_features <- length(X)
    
    if (sum(diff(sapply(X, nrow)))!=0) stop("Number of rows of X should be equal.")
    if((length(Y) != nrow(X[[1]]))) stop("Number of rows of Y and X should be equal.")
    
    
    if (na.rm == TRUE) {
      # NAs are removed here:
      ok = complete.cases(Y,X)
      X <- lapply(X, function(x){x[ok,]})
      Y = Y[ok]
    }
    
    
    n <- length(Y)
    
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    p = length(X)
    
    if (num_features > p) stop("Number of features should not be larger than maximum number of original features.")
    if ((floor(num_features) != num_features) || (num_features <= 0)) stop("Number of features should be a positive integer.")
    
    if (!is.numeric(Y)) stop("currently FOCI does not handle factor Y")
    
    # if (standardize == "scale") {
    #   for (i in 1:p) {
    #     if(length(unique(X[, i])) > 1) {
    #       X[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])
    #     }else{
    #       # NM, May 12; changed to paste0() to remove superfluous blank
    #       stop(paste0("Column ", i, " of X is constant."))
    #     }
    #   }
    # }
    # 
    # if (standardize == "bounded") {
    #   for (i in 1:p) {
    #     if(length(unique(X[, i])) > 1) {
    #       X[,i] <- (X[,i] - min(X[,i])) / (max(X[, i]) - min(X[, i]))
    #     }else{
    #       stop(paste0("Column ", i, " of X is constant."))
    #     }
    #   }
    # }
    
    method=method[1]
    
    if (parPlat[1] == 'none') {
      res <- fofcd_main(Y, X, num_features = num_features,
                        stop = stop, numCores = numCores, method = method, control = control)
      feature <- data.frame(col_remain = res$selectedVar$index, stepT = res$stepT)
      
      core_time <- difftime(Sys.time(), core_time, units = 'sec')
      results <- list(feature = feature, 
                      time = core_time)
      return(results)
    }
    
    # NM, May 12:  many datasets are ordered by one or more columns; to
    # preserve iid-ness, randomize the row order; if we get here, we will
    # be chunking
    nr <- nrow(X[[1]])
    permRowNums <- sample(1:nr,nr,replace=FALSE)
    X <- lapply(X, function(x){x[permRowNums,]})
    Y <- Y[permRowNums]
    
    rowNums <- parallel::splitIndices(length(Y), numCores)
    selectFromChunk <- function(nodeNum) {
      myRows <- rowNums[[nodeNum]]
      sel <- fofcd_main(Y[myRows], lapply(X, function(x){x[myRows,]}), stop = stop,
                        numCores = numCores)$selectedVar$index
    }
    
    if(inherits(parPlat,'cluster')) {
      cls <- parPlat
    }else if(parPlat == 'locThreads') {
      # set up the cluster (in multicore case, it's virtual)
      cls <- parallel::makeCluster(numCores)
    } else stop('invalid parPlat')
    
    # worker nodes load library
    parallel::clusterEvalQ(cls, library(FOCI))
    # ship data to workers
    parallel::clusterExport(cls, c('Y', 'X', 'rowNums', 'selectFromChunk'),
                            envir = environment())
    # drop-in replacement for mclapply
    slc <- parallel::parLapply(cls, seq(1, length(cls)), selectFromChunk)
    
    if (printIntermed) print(slc)
    
    slc <- Reduce(union, slc)
    # May 15, Mona removed FOCI:::
    res <- fofcd_main(Y, X[slc], num_features, stop)
    # must translate indices in reduced system to those of original
    newIdxs <- res$selectedVar$index
    origIdxs <- slc[newIdxs]
    res$selectedVar$index <- origIdxs
    
    res$stepT = res$stepT[1:num_features]
    parallel::stopCluster(cls)
    
    feature <- data.frame(col_remain = res$selectedVar$index, stepT = res$stepT)
    
    core_time <- difftime(Sys.time(), core_time, units = 'sec')
    results <- list(feature = feature, 
                    time = core_time)
    return(results)
  }
}


#####################################################################
basis_expansion <- function (X, method = c('fourier_basis', 'ospline_basis', 'kernel'), 
                             control = list()){
  p = length(X)
  ndimX = sapply(X, ncol)
  idxFD = seq_along(X)[ndimX > 1]
  
  if (is.null(names(X))) {
    names(X) <- paste0('x',1:length(X))
    warning('X lacked column names, has been assigned x1, x2,...')
  }
  namesX <- names(X)
  
  if (method == 'kernel'){
    con = list(t = seq(0,1,len=max(ndimX, na.rm = T)))
    con[(namc <- names(control))] <- control
    t = con$t
    
    if(!is.list(t)){
      tRange = range(t)
      t = lapply(X,function(i){
        if(ncol(i) == length(t)){
          return(t)
        }
        if(ncol(i) != length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i) == 1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t) != p)
        t = lapply(X,function(i)return(t[[1]]))
    }
    
    obs_X <- list()
    approx_X <- list()
    for(iL in idxFD){
      tRange = range(t[[iL]])
      XiL <- X[[iL]]
      obs_X[[iL]] <- matrix(t((apply(XiL, 1, diff)/2 + t(XiL[, 1:(ndimX[iL]-1)]))*
                                diff(t[[iL]]))*(tRange[2]-tRange[1]), ncol = ndimX[[iL]]-1)
      approx_X[[iL]] <- XiL
    }
    names(obs_X) <- namesX
    names(approx_X) <- namesX
  }else
    if (method == 'fourier_basis'){
      stopifnot(require(fda))
      con = list(nbasis = 18, nderiv = 0,
                 t = seq(0, 1, len = max(ndimX, na.rm = T)))
      con[(namc <- names(control))] <- control
      if(is.null(con$samplet)) con$samplet <- con$t
      
      nbasis = con$nbasis
      nderiv = con$nderiv
      t=con$t
      samplet = con$samplet
      
      if(!is.list(t)){
        tRange = range(t)
        t = lapply(X,function(i){
          if(ncol(i) == length(t)){
            return(t)
          }
          if(ncol(i) != length(t) & ncol(i)>1){
            return(seq(tRange[1],tRange[2],len=ncol(i)))
          }
          if(ncol(i) == 1){
            return((tRange[1]+tRange[2])/2)
          }
        })
      }
      if(is.list(t)){
        if(length(t) != p)
          t = lapply(X,function(i)return(t[[1]]))
      }
      
      if(!is.list(samplet)){
        samplet = lapply(X, function(i) return(samplet))
      }
      if(is.list(samplet)){
        if(length(samplet) != p)
          samplet = lapply(X,function(i) return(samplet[[1]]))
      }
      
      if(length(nbasis) != length(nderiv)){
        nbasis = nbasis[1]
        nderiv = nderiv[1]
        warning('length of nbasis is different from length of nderiv')
      }
      if(length(nbasis) != p)
        nbasis = rep(nbasis[1], p)
      if(length(nderiv) != p)
        nderiv = rep(nderiv[1], p)
      
      obs_X <- list()
      approx_X <- list()
      Time_X <- list()
      for(iL in idxFD){
        
        nb = nbasis[iL]
        if(2 * (nb %/% 2) == nb) nb <- nb + 1
        no = nderiv[iL]
        tiL = t[[iL]]
        tRange = range(tiL)
        stiL <- samplet[[iL]]
        
        fbasis <- create.fourier.basis(rangeval = tRange, nbasis = nb)
        X_1 <- t(as.matrix(X[[iL]]))
        fd <- Data2fd(argvals = tiL, 
                      y = X_1, 
                      basisobj= fbasis,
                      lambda=0)
        obs_X[[iL]] <- matrix(t(fd$coefs), ncol = nb)
        approx_X[[iL]] <- X[[iL]]
        Time_X[[iL]] <- tiL
        op <- par(mfrow=c(2,1))
        plot(fbasis, main='bases')
        plot(fd[1], main='fit', ylab = namesX[iL])
        points(tiL, X_1[,1])
        par(op)
        
        if(no){
          for(ideriv in 1:no){
            fd_deriv = deriv.fd(fd, ideriv)
            obs_X[[length(idxFD)*ideriv+iL]] <- matrix(t(fd_deriv$coefs), ncol = nb)
            approx_X[[length(idxFD)*ideriv+iL]] <- matrix(t(eval.fd(stiL, fd_deriv)), 
                                                          ncol = length(stiL))
            Time_X[[length(idxFD)*ideriv+iL]] <- stiL
          }
        }
      }
      namesobs <- c(namesX[idxFD],
                    paste0(namesX[nderiv[idxFD]>=1], rep("_1", times = length(which(nderiv[idxFD]>=1)))),
                    paste0(namesX[nderiv[idxFD]>=2], rep("_2", times = length(which(nderiv[idxFD]>=2))))
      )
      names(obs_X) <- namesobs
      names(approx_X) <- namesobs
      names(Time_X) <- namesobs
    }else
      if (method == 'ospline_basis'){
        stopifnot(require(orthogonalsplinebasis))
        con = list(knots = seq(0, 1, length.out = 10),
                   norder = 4,
                   nderiv = 0,
                   t = seq(0, 1, len = max(ndimX, na.rm = T)))
        con[(namc <- names(control))] <- control
        if(is.null(con$samplet)) con$samplet <- con$t
        knots = con$knots
        norder = con$norder
        nderiv <- con$nderiv
        t=con$t
        samplet = con$samplet
        
        if(!is.list(knots)) knots = lapply(X, function(i)return(knots))
        if(is.list(knots)){
          if(length(knots) != p)
            knots = lapply(X, function(i)return(knots[[1]]))
        }
        
        if(length(norder) != length(nderiv)){
          norder = norder[1]
          nderiv = nderiv[1]
          warning('length of norder is different from length of nderiv')
        }
        if(length(norder) != p)
          norder = rep(norder[1], p)
        if(length(nderiv) != p)
          nderiv = rep(nderiv[1], p)
        
        if(!is.list(t)){
          tRange = range(t)
          t = lapply(X,function(i){
            if(ncol(i) == length(t)){
              return(t)
            }
            if(ncol(i) != length(t) & ncol(i)>1){
              return(seq(tRange[1],tRange[2],len=ncol(i)))
            }
            if(ncol(i) == 1){
              return((tRange[1]+tRange[2])/2)
            }
          })
        }
        if(is.list(t)){
          if(length(t) != p)
            t = lapply(X,function(i)return(t[[1]]))
        }

        
        if(!is.list(samplet)){
          samplet = lapply(X, function(i) return(samplet))
        }
        if(is.list(samplet)){
          if(length(samplet) != p)
            samplet = lapply(X,function(i) return(samplet[[1]]))
        }
        

        
        obs_X <- list()
        approx_X <- list()
        Time_X <- list()
        for(iL in idxFD){
          kn = knots[[iL]]
          no = norder[iL]
          nd = nderiv[iL]
          tiL = t[[iL]]
          stiL <- samplet[[iL]]

          osb <- OrthogonalSplineBasis(knots = expand.knots(kn, order = no))
          csb <- apply(X[[iL]], 1, function(x) fitLS(osb, tiL, as.numeric(x)))
          # approx <- evaluate(osb, tiL) %*% csb
          
          obs_X[[iL]] <- matrix(t(csb), ncol = length(kn)+no-2)
          approx_X[[iL]] <- X[[iL]]
          Time_X[[iL]] <- tiL
          
          op <- par(mfrow=c(2+nd,1))
          plot(osb, main='bases')
          plot(osb, csb[,1], main='fit', ylab = namesX[iL])
          points(tiL, X[[iL]][1,])
          
          if(nd){
            for(ideriv in 1:nd){
              osb_deriv = deriv(osb, ideriv)
              approx <- evaluate(osb_deriv, stiL) %*% csb
              csb_deriv <- apply(t(approx), 1, function(x) fitLS(osb, stiL, as.numeric(x)))
              
              obs_X[[length(idxFD)*ideriv+iL]] <- matrix(t(csb_deriv), ncol = length(kn)+no-2)
              approx_X[[length(idxFD)*ideriv+iL]] <- matrix(t(approx), ncol = length(stiL))
              Time_X[[length(idxFD)*ideriv+iL]] <- stiL
            
              plot(stiL, approx[,1], main='fit', ylab = paste0(namesX[iL], "_", ideriv))
              abline(0,0)
              }
          }
          par(op)
        }
        namesobs <- c(namesX[idxFD],
                      paste0(namesX[nderiv[idxFD]>=1], rep("_1", times = length(which(nderiv[idxFD]>=1)))),
                      paste0(namesX[nderiv[idxFD]>=2], rep("_2", times = length(which(nderiv[idxFD]>=2))))
        )
        names(obs_X) <- namesobs
        names(approx_X) <- namesobs
        names(Time_X) <- namesobs
        
      }else{
        stop("the avaiable method till now is 'raw', fourier basis, and ospline basis.")
      }
  basis_X <- list(obs_X = obs_X, approx_X = approx_X, Time_X = Time_X)
  return(basis_X)
}
