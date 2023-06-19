library(flars)
source("function_FOFCD.R")
source("function_FCCA-DC.R")
library("fda.usc")

iteration = 100

################################################################################
# variable selection 3 - functional linear model
cat("variable selection 3 ", "functional linear model\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100
r <- c(sqrt(0.1), sqrt(0.2), sqrt(0.5))

for (r_idex in 1:length(r)){ 
  result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
  for (iter in 1:iteration){
    cat("\n r^2=", (r[r_idex])^2, ", iteration times ", iter, "\n")
    beta_bar <- c(0.3, 4*(-1)^c(2:D)/(2:D)^2)
    beta <- r[r_idex]*beta_bar/norm(beta_bar, "2")
    phi <- t(cbind(rep(1, times = length(Time)), t(tcrossprod(sqrt(2)*cos(c(1:(D-1))*pi), Time))))
    Beta_t <- crossprod(beta, phi)
    
    X_1 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X1_t <- crossprod(t(X_1), c(1:D)^(-1.1/2)*phi)
    
    X_2 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X2_t <- crossprod(t(X_2), c(1:D)^(-1.1/2)*phi)
    
    X_3 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X3_t <- crossprod(t(X_3), c(1:D)^(-1.1/2)*phi)
    
    X_4 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X4_t <- crossprod(t(X_4), c(1:D)^(-1.1/2)*phi)
    
    f1 <- tcrossprod(X1_t, Beta_t)
    f2 <- tcrossprod(X2_t, Beta_t)
    eps <- rnorm(n, mean = 0, sd = 1)
    y <- f1 + f2 + eps
    
    data <- list(y = matrix(y, ncol = 1), 
                 x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
                 Time = Time)
    
    ## FOFCD
    method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
    fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
    # print(fofcd.fb)
    feature.fb <- fofcd.fb$feature$col_remain
    logic.fb <- (c(1:p) %in% feature.fb)
    if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
    result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
    
    method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                       tmax, tmax),
                                             norder = 2, t = Time)
    fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
    # print(fofcd.osb)
    feature.osb <- fofcd.osb$feature$col_remain
    logic.osb <- (c(1:p) %in% feature.osb)
    if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
    result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
    
    
    # fLARS
    time.fLARS <- Sys.time()
    flars.basis <- flars::flars(data$x, data$y, method='basis', 
                                control = list(t = Time), 
                                max_selection = length(data$x),
                                VarThreshold = 0.9,
                                normalize='norm', lasso=FALSE)
    time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
    feature.fLARS <- flars.basis$index
    # print(feature.fLARS)
    logic.fLARS <- (c(1:p) %in% feature.fLARS)
    if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
    result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
    
    ## VS for glm
    x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
    x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
    x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
    x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
    dat <- data.frame("y"=y)
    ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
    time.glm <- Sys.time()
    res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
    if('try-error' %in% class(res.glm)) {
      feature.glm <- NULL
    }else{
      feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
    }
    time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
    logic.glm <- (c(1:p) %in% feature.glm)
    if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
    result.glm <- rbind(result.glm, c(logic.glm, time.glm))
    rm(ldat)
    
    ## VS for gam
    x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
    x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
    x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
    x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
    dat <- data.frame("y"=y)
    ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
    time.gam <- Sys.time()
    res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
    if('try-error' %in% class(res.gam)) {
      feature.gam <- NULL
    }else{
      feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
    }
    time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
    logic.gam <- (c(1:p) %in% feature.gam)
    if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
    result.gam <- rbind(result.gam, c(logic.gam, time.gam))
    rm(ldat)
  }
  
  colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
    colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
  write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_linear_", 
                                    (r[r_idex])^2*10, "_fofcd_fb.csv"), row.names = F)
  write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_linear_", 
                                     (r[r_idex])^2*10, "_fofcd_osb.csv"), row.names = F)
  write.csv(result.flars.basis, paste0("result_variable_selection_simu_linear_", 
                                       (r[r_idex])^2*10, "_flars_basis.csv"), row.names = F)
  write.csv(result.glm, paste0("result_variable_selection_simu_linear_", 
                               (r[r_idex])^2*10, "_glm.csv"), row.names = F)
  write.csv(result.gam, paste0("result_variable_selection_simu_linear_", 
                               (r[r_idex])^2*10, "_gam.csv"), row.names = F)						   
}


################################################################################
# variable selection 4 - functional linear model with quadratic items
cat("variable selection 4 ", "functional linear model with quadratic items\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100
r <- c(0.5, 1, 2)

for (r_idex in 1:length(r)){ 
  result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
  for (iter in 1:iteration){
    cat("\n r=", r[r_idex], ", iteration times ", iter, "\n")
    phi <- t(cbind(rep(1, times = length(Time)), t(tcrossprod(sqrt(2)*cos(c(1:(D-1))*pi), Time))))
    
    X_1 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X1_t <- crossprod(t(X_1), c(1:D)^(-1.1/2)*phi)
    
    X_2 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X2_t <- crossprod(t(X_2), c(1:D)^(-1.1/2)*phi)
    
    X_3 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X3_t <- crossprod(t(X_3), c(1:D)^(-1.1/2)*phi)
    
    X_4 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X4_t <- crossprod(t(X_4), c(1:D)^(-1.1/2)*phi)
    
    f1 <- r[r_idex]*apply(X1_t, 1, crossprod)
    f2 <- r[r_idex]*apply(X2_t, 1, crossprod)
    eps <- rnorm(n, mean = 0, sd = 1)
    y <- f1 + f2 + eps
    
    data <- list(y = matrix(y, ncol = 1), 
                 x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
                 Time = Time)
    
    ## FOFCD
    method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
    fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
    # print(fofcd.fb)
    feature.fb <- fofcd.fb$feature$col_remain
    logic.fb <- (c(1:p) %in% feature.fb)
    if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
    result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
    
    method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                       tmax, tmax),
                                             norder = 2, t = Time)
    fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
    # print(fofcd.osb)
    feature.osb <- fofcd.osb$feature$col_remain
    logic.osb <- (c(1:p) %in% feature.osb)
    if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
    result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
    
    
    # fLARS
    time.fLARS <- Sys.time()
    flars.basis <- flars::flars(data$x, data$y, method='basis', 
                                control = list(t = Time), 
                                max_selection = length(data$x),
                                VarThreshold = 0.9,
                                normalize='norm', lasso=FALSE)
    time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
    feature.fLARS <- flars.basis$index
    # print(feature.fLARS)
    logic.fLARS <- (c(1:p) %in% feature.fLARS)
    if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
    result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
    
    ## VS for glm
    x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
    x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
    x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
    x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
    dat <- data.frame("y"=y)
    ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
    time.glm <- Sys.time()
    res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
    if('try-error' %in% class(res.glm)) {
      feature.glm <- NULL
    }else{
      feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
    }
    time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
    logic.glm <- (c(1:p) %in% feature.glm)
    if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
    result.glm <- rbind(result.glm, c(logic.glm, time.glm))
    rm(ldat)
    
    ## VS for gam
    x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
    x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
    x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
    x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
    dat <- data.frame("y"=y)
    ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
    time.gam <- Sys.time()
    res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
    if('try-error' %in% class(res.gam)) {
      feature.gam <- NULL
    }else{
      feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
    }
    time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
    logic.gam <- (c(1:p) %in% feature.gam)
    if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
    result.gam <- rbind(result.gam, c(logic.gam, time.gam))
    rm(ldat)
  }
  colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
    colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
  write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_quadratic_", 
                                    r[r_idex]*10, "_fofcd_fb.csv"), row.names = F)
  write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_quadratic_", 
                                     r[r_idex]*10, "_fofcd_osb.csv"), row.names = F)
  write.csv(result.flars.basis, paste0("result_variable_selection_simu_quadratic_", 
                                       r[r_idex]*10, "_flars_basis.csv"), row.names = F)
  write.csv(result.glm, paste0("result_variable_selection_simu_quadratic_", 
                               r[r_idex]*10, "_glm.csv"), row.names = F)	
  write.csv(result.gam, paste0("result_variable_selection_simu_quadratic_", 
                               r[r_idex]*10, "_gam.csv"), row.names = F)											
}


################################################################################
# variable selection 5.1 - functional additive model with quadratic items
cat("variable selection 5.1 ", "functional additive model with quadratic items\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100
r <- c(0.1, 0.2, 0.3)

for (r_idex in 1:length(r)){ 
  result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
  result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
  for (iter in 1:iteration){
    cat("\n r=", r[r_idex], ", iteration times ", iter, "\n")
    phi <- t(cbind(rep(1, times = length(Time)), t(tcrossprod(sqrt(2)*cos(c(1:(D-1))*pi), Time))))
    
    X_1 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X1_t <- crossprod(t(X_1), c(1:D)^(-1.1/2)*phi)
    
    X_2 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X2_t <- crossprod(t(X_2), c(1:D)^(-1.1/2)*phi)
    
    X_3 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X3_t <- crossprod(t(X_3), c(1:D)^(-1.1/2)*phi)
    
    X_4 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
    X4_t <- crossprod(t(X_4), c(1:D)^(-1.1/2)*phi)
    
    f1 <- apply(X1_t, 1, crossprod)
    f2 <- apply(X2_t, 1, crossprod)
    eps <- rnorm(n, mean = 0, sd = 1)
    y <- log(f1) + r[r_idex]*(log(f2))^2 + eps
    
    data <- list(y = matrix(y, ncol = 1), 
                 x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
                 Time = Time)
    
    ## FOFCD
    method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
    fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
    # print(fofcd.fb)
    feature.fb <- fofcd.fb$feature$col_remain
    logic.fb <- (c(1:p) %in% feature.fb)
    if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
    result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
    
    method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                       tmax, tmax),
                                             norder = 2, t = Time)
    fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
    # print(fofcd.osb)
    feature.osb <- fofcd.osb$feature$col_remain
    logic.osb <- (c(1:p) %in% feature.osb)
    if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
    result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
    
    
    # fLARS
    time.fLARS <- Sys.time()
    flars.basis <- flars::flars(data$x, data$y, method='basis', 
                                control = list(t = Time), 
                                max_selection = length(data$x),
                                VarThreshold = 0.9,
                                normalize='norm', lasso=FALSE)
    time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
    feature.fLARS <- flars.basis$index
    # print(feature.fLARS)
    logic.fLARS <- (c(1:p) %in% feature.fLARS)
    if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
    result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
    
    ## VS for glm
    x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
    x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
    x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
    x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
    dat <- data.frame("y"=y)
    ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
    time.glm <- Sys.time()
    res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
    if('try-error' %in% class(res.glm)) {
      feature.glm <- NULL
    }else{
      feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
    }
    time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
    logic.glm <- (c(1:p) %in% feature.glm)
    if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
    result.glm <- rbind(result.glm, c(logic.glm, time.glm))
    rm(ldat)
    
    ## VS for gam
    x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
    x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
    x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
    x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
    dat <- data.frame("y"=y)
    ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
    time.gam <- Sys.time()
    res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
    if('try-error' %in% class(res.gam)) {
      feature.gam <- NULL
    }else{
      feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
    }
    time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
    logic.gam <- (c(1:p) %in% feature.gam)
    if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
    result.gam <- rbind(result.gam, c(logic.gam, time.gam))
    rm(ldat)
  }
  colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
    colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
  write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_additive1_", 
                                    r[r_idex]*10, "_fofcd_fb.csv"), row.names = F)
  write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_additive1_", 
                                     r[r_idex]*10, "_fofcd_osb.csv"), row.names = F)
  write.csv(result.flars.basis, paste0("result_variable_selection_simu_additive1_", 
                                       r[r_idex]*10, "_flars_basis.csv"), row.names = F)
  write.csv(result.glm, paste0("result_variable_selection_simu_additive1_", 
                               r[r_idex]*10, "_glm.csv"), row.names = F)
  write.csv(result.gam, paste0("result_variable_selection_simu_additive1_", 
                               r[r_idex]*10, "_gam.csv"), row.names = F)
}

################################################################################
# variable selection 5.2 - functional additive model with quadratic items
cat("variable selection 5.2 ", "functional additive model with quadratic items\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100

result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
for (iter in 1:iteration){
  cat("\n iteration times ", iter, "\n")
  phi <- t(cbind(rep(1, times = length(Time)), t(tcrossprod(sqrt(2)*cos(c(1:(D-1))*pi), Time))))
  
  X_1 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X1_t <- crossprod(t(X_1), c(1:D)^(-1.1/2)*phi)
  
  X_2 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X2_t <- crossprod(t(X_2), c(1:D)^(-1.1/2)*phi)
  
  X_3 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X3_t <- crossprod(t(X_3), c(1:D)^(-1.1/2)*phi)
  
  X_4 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X4_t <- crossprod(t(X_4), c(1:D)^(-1.1/2)*phi)
  
  f1 <- apply(X1_t, 1, crossprod)
  f2 <- apply(X2_t, 1, crossprod)
  eps <- rnorm(n, mean = 0, sd = 1)
  y <- log(f1) + log(f2) + eps
  
  data <- list(y = matrix(y, ncol = 1), 
               x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
               Time = Time)
  
  ## FOFCD
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.fb)
  feature.fb <- fofcd.fb$feature$col_remain
  logic.fb <- (c(1:p) %in% feature.fb)
  if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
  result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.osb)
  feature.osb <- fofcd.osb$feature$col_remain
  logic.osb <- (c(1:p) %in% feature.osb)
  if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
  result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
  
  
  # fLARS
  time.fLARS <- Sys.time()
  flars.basis <- flars::flars(data$x, data$y, method='basis', 
                              control = list(t = Time), 
                              max_selection = length(data$x),
                              VarThreshold = 0.9,
                              normalize='norm', lasso=FALSE)
  time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
  feature.fLARS <- flars.basis$index
  # print(feature.fLARS)
  logic.fLARS <- (c(1:p) %in% feature.fLARS)
  if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
  result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
  
  ## VS for glm
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.glm <- Sys.time()
  res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.glm)) {
    feature.glm <- NULL
  }else{
    feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
  }
  time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
  logic.glm <- (c(1:p) %in% feature.glm)
  if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
  result.glm <- rbind(result.glm, c(logic.glm, time.glm))
  rm(ldat)
  
  ## VS for gam
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.gam <- Sys.time()
  res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.gam)) {
    feature.gam <- NULL
  }else{
    feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
  }
  time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
  logic.gam <- (c(1:p) %in% feature.gam)
  if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
  result.gam <- rbind(result.gam, c(logic.gam, time.gam))
  rm(ldat)
}
colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
  colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_additive2_fofcd_fb.csv"), row.names = F)
write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_additive2_fofcd_osb.csv"), row.names = F)
write.csv(result.flars.basis, paste0("result_variable_selection_simu_additive2_flars_basis.csv"), row.names = F)
write.csv(result.glm, paste0("result_variable_selection_simu_additive2_glm.csv"), row.names = F)
write.csv(result.gam, paste0("result_variable_selection_simu_additive2_gam.csv"), row.names = F)


################################################################################
# variable selection 5.3 - functional additive model with quadratic items
cat("variable selection 5.3 ", "functional additive model with quadratic items\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100

result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
for (iter in 1:iteration){
  cat("\n iteration times ", iter, "\n")
  phi <- t(cbind(rep(1, times = length(Time)), t(tcrossprod(sqrt(2)*cos(c(1:(D-1))*pi), Time))))
  
  X_1 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X1_t <- crossprod(t(X_1), c(1:D)^(-1.1/2)*phi)
  
  X_2 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X2_t <- crossprod(t(X_2), c(1:D)^(-1.1/2)*phi)
  
  X_3 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X3_t <- crossprod(t(X_3), c(1:D)^(-1.1/2)*phi)
  
  X_4 <- matrix(rnorm(n*D, mean = 0, sd = 1), nrow = n)
  X4_t <- crossprod(t(X_4), c(1:D)^(-1.1/2)*phi)
  
  f1 <- apply(X1_t, 1, crossprod)
  f2 <- apply(X2_t, 1, crossprod)
  eps <- rnorm(n, mean = 0, sd = 1)
  y <- log(f1 + f2) + eps
  
  data <- list(y = matrix(y, ncol = 1), 
               x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
               Time = Time)
  
  ## FOFCD
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.fb)
  feature.fb <- fofcd.fb$feature$col_remain
  logic.fb <- (c(1:p) %in% feature.fb)
  if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
  result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.osb)
  feature.osb <- fofcd.osb$feature$col_remain
  logic.osb <- (c(1:p) %in% feature.osb)
  if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
  result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
  
  
  # fLARS
  time.fLARS <- Sys.time()
  flars.basis <- flars::flars(data$x, data$y, method='basis', 
                              control = list(t = Time), 
                              max_selection = length(data$x),
                              VarThreshold = 0.9,
                              normalize='norm', lasso=FALSE)
  time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
  feature.fLARS <- flars.basis$index
  # print(feature.fLARS)
  logic.fLARS <- (c(1:p) %in% feature.fLARS)
  if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
  result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
  
  ## VS for glm
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.glm <- Sys.time()
  res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.glm)) {
    feature.glm <- NULL
  }else{
    feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
  }
  time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
  logic.glm <- (c(1:p) %in% feature.glm)
  if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
  result.glm <- rbind(result.glm, c(logic.glm, time.glm))
  rm(ldat)
  
  ## VS for gam
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.gam <- Sys.time()
  res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.gam)) {
    feature.gam <- NULL
  }else{
    feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
  }
  time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
  logic.gam <- (c(1:p) %in% feature.gam)
  if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
  result.gam <- rbind(result.gam, c(logic.gam, time.gam))
  rm(ldat)
}
colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
  colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_additive3_fofcd_fb.csv"), row.names = F)
write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_additive3_fofcd_osb.csv"), row.names = F)
write.csv(result.flars.basis, paste0("result_variable_selection_simu_additive3_flars_basis.csv"), row.names = F)
write.csv(result.glm, paste0("result_variable_selection_simu_additive3_glm.csv"), row.names = F)
write.csv(result.gam, paste0("result_variable_selection_simu_additive3_gam.csv"), row.names = F)



################################################################################
# variable selection 6.1 - functional linear model with interactions
cat("variable selection 6.1 ", "functional linear model with interactions\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100
result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
for (iter in 1:iteration){
  cat("\n iteration times ", iter, "\n")
  phi <- t(cbind(sqrt(2)*cos(pi*Time), sqrt(2)*sin(pi*Time),
                 sqrt(2)*cos(2*pi*Time), sqrt(2)*sin(2*pi*Time),
                 sqrt(2)*cos(3*pi*Time), sqrt(2)*sin(3*pi*Time)))
  tau <- matrix(c(4, 0.5, 0, 0, 4, 1, 0.5, 0, 4, 2, 1, 0, 4, 1, 0, 0), nrow = 4, byrow = TRUE)
  
  xi <- array(0, dim = c(n, 4, 4))
  for (ip in 1:nrow(tau)){
    for (iq in 1:ncol(tau)){
      xi[, ip, iq] <- rnorm(n, mean = 0, sd = sqrt(tau[ip, iq]))
    }
  }
  X1_t <- crossprod(t(xi[,1,1:2]), phi[1:2,])
  X2_t <- crossprod(t(xi[,2,1:3]), phi[1:3,])
  X3_t <- crossprod(t(xi[,3,1:3]), phi[c(1,4,6),])
  X4_t <- crossprod(t(xi[,4,1:2]), phi[5:6,])
  
  f1 <- xi[,1,2]*xi[,2,1]
  eps <- rnorm(n, mean = 0, sd = 1)
  y <- f1 + eps
  
  data <- list(y = matrix(y, ncol = 1), 
               x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
               Time = Time)
  
  ## FOFCD
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.fb)
  feature.fb <- fofcd.fb$feature$col_remain
  logic.fb <- (c(1:p) %in% feature.fb)
  if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
  result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.osb)
  feature.osb <- fofcd.osb$feature$col_remain
  logic.osb <- (c(1:p) %in% feature.osb)
  if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
  result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
  
  
  # fLARS
  time.fLARS <- Sys.time()
  flars.basis <- flars::flars(data$x, data$y, method='basis', 
                              control = list(t = Time), 
                              max_selection = length(data$x),
                              VarThreshold = 0.9,
                              normalize='norm', lasso=FALSE)
  time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
  feature.fLARS <- flars.basis$index
  # print(feature.fLARS)
  logic.fLARS <- (c(1:p) %in% feature.fLARS)
  if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
  result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
  
  ## VS for glm
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.glm <- Sys.time()
  res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.glm)) {
    feature.glm <- NULL
  }else{
    feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
  }
  time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
  logic.glm <- (c(1:p) %in% feature.glm)
  if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
  result.glm <- rbind(result.glm, c(logic.glm, time.glm))
  rm(ldat)
  
  ## VS for gam
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.gam <- Sys.time()
  res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.gam)) {
    feature.gam <- NULL
  }else{
    feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
  }
  time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
  logic.gam <- (c(1:p) %in% feature.gam)
  if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
  result.gam <- rbind(result.gam, c(logic.gam, time.gam))
  rm(ldat)
}
colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
  colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_interaction1_fofcd_fb.csv"), row.names = F)
write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_interaction1_fofcd_osb.csv"), row.names = F)
write.csv(result.flars.basis, paste0("result_variable_selection_simu_interaction1_flars_basis.csv"), row.names = F)
write.csv(result.glm, paste0("result_variable_selection_simu_interaction1_glm.csv"), row.names = F)
write.csv(result.gam, paste0("result_variable_selection_simu_interaction1_gam.csv"), row.names = F)


################################################################################
# variable selection 6.2 - functional linear model with interactions
cat("variable selection 6.2 ", "functional linear model with interactions\n")
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)
p <- 4
D <- 100
result.fofcd.fb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.fofcd.osb <- data.frame(matrix(ncol = p+1, nrow = 0))
result.flars.basis <- data.frame(matrix(ncol = p+1, nrow = 0))
result.glm <- data.frame(matrix(ncol = p+1, nrow = 0))
result.gam <- data.frame(matrix(ncol = p+1, nrow = 0))
for (iter in 1:iteration){
  cat("\n iteration times ", iter, "\n")
  phi <- t(cbind(sqrt(2)*cos(pi*Time), sqrt(2)*sin(pi*Time),
                 sqrt(2)*cos(2*pi*Time), sqrt(2)*sin(2*pi*Time),
                 sqrt(2)*cos(3*pi*Time), sqrt(2)*sin(3*pi*Time)))
  tau <- matrix(c(4, 0.5, 0, 0, 4, 1, 0.5, 0, 4, 2, 1, 0, 4, 1, 0, 0), nrow = 4, byrow = TRUE)
  
  xi <- array(0, dim = c(n, 4, 4))
  for (ip in 1:nrow(tau)){
    for (iq in 1:ncol(tau)){
      xi[, ip, iq] <- rnorm(n, mean = 0, sd = sqrt(tau[ip, iq]))
    }
  }
  X1_t <- crossprod(t(xi[,1,1:2]), phi[1:2,])
  X2_t <- crossprod(t(xi[,2,1:3]), phi[1:3,])
  X3_t <- crossprod(t(xi[,3,1:3]), phi[c(1,4,6),])
  X4_t <- crossprod(t(xi[,4,1:2]), phi[5:6,])
  
  f1 <- xi[,1,2]*xi[,2,1] + xi[,1,1]*xi[,2,3] # xi[,1,1]*xi[,2,2]/2
  eps <- rnorm(n, mean = 0, sd = 1)
  y <- f1 + eps
  
  data <- list(y = matrix(y, ncol = 1), 
               x = list(x1 = X1_t, x2 = X2_t, x3 = X3_t, x4 = X4_t), 
               Time = Time)
  
  ## FOFCD
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.fb)
  feature.fb <- fofcd.fb$feature$col_remain
  logic.fb <- (c(1:p) %in% feature.fb)
  if (!is.null(feature.fb)) logic.fb[which(logic.fb==1)] <- order(feature.fb)
  result.fofcd.fb <- rbind(result.fofcd.fb, c(logic.fb, fofcd.fb$time))
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  fofcd.osb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.osb)
  feature.osb <- fofcd.osb$feature$col_remain
  logic.osb <- (c(1:p) %in% feature.osb)
  if (!is.null(feature.osb)) logic.osb[which(logic.osb==1)] <- order(feature.osb)
  result.fofcd.osb <- rbind(result.fofcd.osb, c(logic.osb, fofcd.osb$time))
  
  
  # fLARS
  time.fLARS <- Sys.time()
  flars.basis <- flars::flars(data$x, data$y, method='basis', 
                              control = list(t = Time), 
                              max_selection = length(data$x),
                              VarThreshold = 0.9,
                              normalize='norm', lasso=FALSE)
  time.fLARS <- difftime(Sys.time(), time.fLARS, units = 'sec')
  feature.fLARS <- flars.basis$index
  # print(feature.fLARS)
  logic.fLARS <- (c(1:p) %in% feature.fLARS)
  if (!is.null(feature.fLARS)) logic.fLARS[which(logic.fLARS==1)] <- order(feature.fLARS)
  result.flars.basis <- rbind(result.flars.basis, c(logic.fLARS, time.fLARS))
  
  
  ## VS for glm
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.glm <- Sys.time()
  res.glm <- try(fregre.glm.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.glm)) {
    feature.glm <- NULL
  }else{
    feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
  }
  time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
  logic.glm <- (c(1:p) %in% feature.glm)
  if (!is.null(feature.glm)) logic.glm[which(logic.glm==1)] <- order(feature.glm)
  result.glm <- rbind(result.glm, c(logic.glm, time.glm))
  rm(ldat)
  
  ## VS for gam
  x1 <- fdata(X1_t, argvals = Time, rangeval = TRange)
  x2 <- fdata(X2_t, argvals = Time, rangeval = TRange)
  x3 <- fdata(X3_t, argvals = Time, rangeval = TRange)
  x4 <- fdata(X4_t, argvals = Time, rangeval = TRange)
  dat <- data.frame("y"=y)
  ldat <- ldata(df=dat, "x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4)
  time.gam <- Sys.time()
  res.gam <- try(fregre.gsam.vs(data=ldat,y="y",numbasis.opt=T, dcor.min=0.01))
  if('try-error' %in% class(res.gam)) {
    feature.gam <- NULL
  }else{
    feature.gam <- as.numeric(stringr::str_remove(res.gam$ipredictors, pattern = "x"))
  }
  time.gam <- difftime(Sys.time(), time.gam, units = 'sec')
  logic.gam <- (c(1:p) %in% feature.gam)
  if (!is.null(feature.gam)) logic.gam[which(logic.gam==1)] <- order(feature.gam)
  result.gam <- rbind(result.gam, c(logic.gam, time.gam))
  rm(ldat)
}
colnames(result.fofcd.fb) <- colnames(result.fofcd.osb) <- colnames(result.glm) <- 
  colnames(result.gam) <- colnames(result.flars.basis) <- c(1:p, "time")
write.csv(result.fofcd.fb, paste0("result_variable_selection_simu_interaction2_fofcd_fb.csv"), row.names = F)
write.csv(result.fofcd.osb, paste0("result_variable_selection_simu_interaction2_fofcd_osb.csv"), row.names = F)
write.csv(result.flars.basis, paste0("result_variable_selection_simu_interaction2_flars_basis.csv"), row.names = F)
write.csv(result.glm, paste0("result_variable_selection_simu_interaction2_glm.csv"), row.names = F)
write.csv(result.gam, paste0("result_variable_selection_simu_interaction2_gam.csv"), row.names = F)

