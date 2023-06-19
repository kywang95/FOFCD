##################################################
#                                                #
##           Functional Weight Plot              #
## Beijing Multi-Site Air-Quality Data Data Set  #
#                                                #
##################################################

# Libraries
library(fda)
library(mice)
library(flars)
library(fda.usc)
library(stringr)
source("function_FOFCD.R")
source("function_FCCA-DC.R")

# Loading data
path <- ("../函数型数据特征选择方法-数据/PRSA_Data_20130301-20170228/")
locations <- c("Guanyuan")
for (iloc in 1:length(locations)){
  location <- locations[iloc]
  data <- read.csv(paste0(path, "PRSA_Data_", location, "_20130301-20170228.csv"), header = T)
  data.2015 <- data[which(data$year=="2015"), ]
  
  n <- nrow(data.2015)/24
  # define the time points on which the functional predictor is observed. 
  timepts <- unique(data.2015$hour)/24
  ntime <- length(timepts)
  Trange <- range(timepts)
  
  repitition <- 100
  result.MSEP <- matrix(nrow = repitition, ncol = 5)
  result.feature <- matrix(nrow = 8*5, ncol = repitition)
  for (irep in 1:repitition){
    cat("repitition ", irep, "\n")
    
    obs <- data.2015[, -c(1:5, ncol(data.2015))]
    y <- matrix(0, nrow = n, ncol = 1)
    scalar_full <- matrix(0, nrow = n, ncol = 1)
    PM10 <- SO2 <- NO2 <- CO <- O3 <- TEMP <- PRES <- DEWP <- RAIN <- 
      wd <- WSPM <- matrix(0, nrow = n, ncol = ntime)
    for(iobs in 1:n){
      # Obtaining response
      y[iobs, ] <- log(mean(obs[(iobs-1)*24+c(1:24), 1], na.rm = T))
      # y[iobs, ] <- log10(mean(obs[(iobs-1)*24+c(1:24), 1], na.rm = T))
      scalar_full[iobs, ] <- sum(obs[(iobs-1)*24+c(1:24), 1], na.rm = T)
      PM10[iobs, ] <- obs[(iobs-1)*24+c(1:24), 2]
      SO2[iobs, ] <- obs[(iobs-1)*24+c(1:24), 3]
      NO2[iobs, ] <- obs[(iobs-1)*24+c(1:24), 4]
      CO[iobs, ] <- obs[(iobs-1)*24+c(1:24), 5]
      O3[iobs, ] <- obs[(iobs-1)*24+c(1:24), 6]
      TEMP[iobs, ] <- obs[(iobs-1)*24+c(1:24), 7]
      PRES[iobs, ] <- obs[(iobs-1)*24+c(1:24), 8]
      DEWP[iobs, ] <- obs[(iobs-1)*24+c(1:24), 9]
      RAIN[iobs, ] <- obs[(iobs-1)*24+c(1:24), 10]
      wd[iobs, ] <- obs[(iobs-1)*24+c(1:24), 11]
      WSPM[iobs, ] <- obs[(iobs-1)*24+c(1:24), 12]
    }
    
    PM10 <- complete(mice(PM10, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    SO2 <- complete(mice(SO2, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    NO2 <- complete(mice(NO2, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    CO <- complete(mice(CO, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    O3 <- complete(mice(O3, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    TEMP <- complete(mice(TEMP, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    PRES <- complete(mice(PRES, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    DEWP <- complete(mice(DEWP, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    RAIN <- complete(mice(RAIN, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    WSPM <- complete(mice(WSPM, m=1, maxit=50, method = 'norm.predict', printFlag = F), action = 1)
    
    data <- list(y=y, 
                 PM10=as.matrix(PM10), SO2=as.matrix(SO2), 
                 NO2=as.matrix(NO2), CO=as.matrix(CO),
                 O3=as.matrix(O3), TEMP=as.matrix(TEMP), 
                 # PRES=as.matrix(PRES), RAIN=as.matrix(RAIN), 
                 DEWP=as.matrix(DEWP), WSPM=as.matrix(WSPM))
    
    #############################################################################
    ## divide the dataset into train set and test set
    #############################################################################
    
    ntrain <- 300
    Train <- lapply(data, function(x) x[1:ntrain,])
    # Train$timepts <- timepts 
    ntest <- dim(data[[2]])[1] - ntrain
    Test <- lapply(data, function(x) x[-c(1:ntrain),])
    # Test$timepts <- timepts
    
    Train_x <- Train[-1]
    Train_y <- matrix(Train[[1]], ncol = 1)
    Test_x <- Test[-1]
    Test_y <- matrix(Test[[1]], ncol = 1)
    colnames(Train_y) <- colnames(Test_y) <- names(data)[1]
    
    ndimX = lapply(Train_x, ncol)
    names(Train_x) <- names(data)[-1]
    names(Test_x) <- names(data)[-1]
    
    # define the fourier basis 
    nbasis = 15
    temp_basis = create.fourier.basis(Trange, nbasis)
    # convert the functional predictor in training data into fda objects
    p <- length(Train_x)
    train_n <- nrow(Train_x[[1]])
    test_n <- nrow(Test_x[[1]])
    train_fda = array(dim = c(nbasis, train_n, p))
    test_fda = array(dim = c(nbasis, test_n, p))
    train_data <- list()
    test_data <- list()
    for(ip in 1:p){
      train_x_data =  Data2fd(timepts, t(Train_x[[ip]]), temp_basis)
      train_data[[ip]] <- t(eval.fd(timepts, train_x_data))
      func_cov = train_x_data$coefs
      train_fda[,,ip] <- func_cov
      test_x_data =  Data2fd(timepts, t(Test_x[[ip]]), temp_basis)
      test_data[[ip]] <- t(eval.fd(timepts, test_x_data))
      func_cov = test_x_data$coefs
      test_fda[,,ip] <- func_cov
    }
    names(train_data) <- paste0("X", c(1:p))
    names(test_data) <- paste0("X", c(1:p))
    
    #############################################################################################
    ## variable selection
    #############################################################################################
    
    ## FOFCD
    method = 'fourier_basis'; control = list(nbasis = 2, norder = 2, t = timepts)
    fofcd.fb <- fofcd(Train_y, train_data, numCores = 1, method = method, control = control)
    # print(fofcd.fb)
    feature.fb <- fofcd.fb$feature$col_remain
    print(feature.fb)
    
    tmax = max(timepts)
    tmin = min(timepts)
    method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                       tmax, tmax),
                                             norder = 2, t = timepts)
    fofcd.osb <- fofcd(Train_y, train_data, numCores = 1, method = method, control = control)
    # print(fofcd.osb)
    feature.osb <- fofcd.osb$feature$col_remain
    print(feature.osb)
    
    ## fLARS
    flars.basis <- flars::flars(train_data, Train_y, method='basis', 
                                control = list(t = timepts), 
                                max_selection = p,
                                normalize='norm'
                                # , lasso=TRUE
    )
    feature.fLARS <- flars.basis$index
    print(feature.fLARS)
    
    #############################################################################
    ## Prediction using fregre.gsam
    #############################################################################
    library(fda.usc)
    x1_train <- fdata(train_data[[1]], argvals = timepts, rangeval = Trange)
    x2_train <- fdata(train_data[[2]], argvals = timepts, rangeval = Trange)
    x3_train <- fdata(train_data[[3]], argvals = timepts, rangeval = Trange)
    x4_train <- fdata(train_data[[4]], argvals = timepts, rangeval = Trange)
    x5_train <- fdata(train_data[[5]], argvals = timepts, rangeval = Trange)
    x6_train <- fdata(train_data[[6]], argvals = timepts, rangeval = Trange)
    x7_train <- fdata(train_data[[7]], argvals = timepts, rangeval = Trange)
    x8_train <- fdata(train_data[[8]], argvals = timepts, rangeval = Trange)
    dat <- data.frame("y"=Train_y)
    ldat <- ldata(df=dat, "x1" = x1_train, "x2" = x2_train, "x3" = x3_train, "x4" = x4_train,
                  "x5" = x5_train, "x6" = x6_train, "x7" = x7_train, "x8" = x8_train)
    
    x1_test <- fdata(test_data[[1]], argvals = timepts, rangeval = Trange)
    x2_test <- fdata(test_data[[2]], argvals = timepts, rangeval = Trange)
    x3_test <- fdata(test_data[[3]], argvals = timepts, rangeval = Trange)
    x4_test <- fdata(test_data[[4]], argvals = timepts, rangeval = Trange)
    x5_test <- fdata(test_data[[5]], argvals = timepts, rangeval = Trange)
    x6_test <- fdata(test_data[[6]], argvals = timepts, rangeval = Trange)
    x7_test <- fdata(test_data[[7]], argvals = timepts, rangeval = Trange)
    x8_test <- fdata(test_data[[8]], argvals = timepts, rangeval = Trange)
    dat <- data.frame("y"=Test_y)
    newldat <- ldata(df=dat, "x1" = x1_test, "x2" = x2_test, "x3" = x3_test, "x4" = x4_test,
                     "x5" = x5_test, "x6" = x6_test, "x7" = x7_test, "x8" = x8_test)
    
    #############################################################################
    ## VS for GLM
    #############################################################################
    time.glm <- Sys.time()
    res.glm <- try(fregre.glm.vs(data=ldat, y="y", numbasis.opt=F
                                  , dcor.min = 0.01
                                  # , dcor.min=1e-8, alpha = 1 ,include = c("x1", "x3", "x6", "x7", "x8")
    ))
    if('try-error' %in% class(res.glm)) {
      feature.glm <- NULL
    }else{
      feature.glm <- as.numeric(stringr::str_remove(res.glm$ipredictors, pattern = "x"))
    }
    time.glm <- difftime(Sys.time(), time.glm, units = 'sec')
    pred.glm.train <- predict(res.glm, ldat)
    print(feature.glm)
    time.glm
    # print(mean((matrix(pred.glm.train, ncol = 1) - Train_y)^2, na.rm = T))
    
    # predict
    pred.glm.test <- predict(res.glm, newldat)
    MSEP.glm <- mean((matrix(pred.glm.test, ncol = 1) - Test_y)^2, na.rm = T)
    
    
    #############################################################################
    ## prediction using features selected by fofcd
    #############################################################################
    features.selected <- paste0("x", feature.fb)
    res.fofcd.fb <- try(fregre.glm.vs(data=ldat, y="y", numbasis.opt=F
                                       , dcor.min=0, alpha = 1
                                       , include = features.selected
                                       # , AIC = FALSE 
    ))
    if('try-error' %in% class(res.fofcd.fb)) {
      feature.fofcd.fb <- NULL
    }else{
      feature.fofcd.fb <- as.numeric(stringr::str_remove(res.fofcd.fb$ipredictors, pattern = "x"))
    }
    pred.fofcd.fb.train <- predict(res.fofcd.fb, ldat)
    feature.fofcd.fb
    # print(mean((matrix(pred.fofcd.fb.train, ncol = 1) - Train_y)^2, na.rm = T))
    
    # predict
    pred.fofcd.fb.test <- predict(res.fofcd.fb, newldat)
    MSEP.fofcd.fb <- mean((matrix(pred.fofcd.fb.test, ncol = 1) - Test_y)^2, na.rm = T)
    
    #############################################################################
    ## prediction using features selected by FOFCD
    #############################################################################
    features.selected <- paste0("x", feature.osb)
    res.fofcd.osb <- try(fregre.glm.vs(data=ldat, y="y", numbasis.opt=F
                                        , dcor.min=0, alpha = 1
                                        , include = features.selected
                                        # , AIC = FALSE
    ))
    if('try-error' %in% class(res.fofcd.osb)) {
      feature.fofcd.osb <- NULL
    }else{
      feature.fofcd.osb <- as.numeric(stringr::str_remove(res.fofcd.osb$ipredictors, pattern = "x"))
    }
    pred.fofcd.osb.train <- predict(res.fofcd.osb, ldat)
    feature.fofcd.osb
    # print(mean((matrix(pred.fofcd.osb.train, ncol = 1) - Train_y)^2, na.rm = T))
    
    # predict
    pred.fofcd.osb.test <- predict(res.fofcd.osb, newldat)
    MSEP.fofcd.osb <- mean((matrix(pred.fofcd.osb.test, ncol = 1) - Test_y)^2, na.rm = T)
    
    #############################################################################
    ## prediction using features selected by fLARS
    #############################################################################
    features.selected <- paste0("x", feature.fLARS)
    res.fLARS <- try(
      fregre.glm.vs(data=ldat, y="y", numbasis.opt=F
                     , dcor.min=0, alpha = 1
                     , include = features.selected
                     # , AIC = FALSE 
                     # , trace = T
      )
    )
    if('try-error' %in% class(res.fLARS)) {
      feature.flars <- NULL
      res.fLARS <- fregre.glm.vs(data=ldat, y="y", numbasis.opt=F
                      , dcor.min=0, alpha = 1
                      , include = features.selected
                      , AIC = FALSE
                      # , trace = T
        )
    }
    feature.flars <- as.numeric(stringr::str_remove(res.fLARS$ipredictors, pattern = "x"))
    # pred.fLARS.train <- predict(res.fLARS, ldat)
    feature.flars
    # print(mean((matrix(pred.fLARS.train, ncol = 1) - Train_y)^2, na.rm = T))
    
    # predict
    if (length(feature.flars) == 0)  {
      pred.fLARS.test <- rep(mean(Test_y), times = length(Test_y))
    }else{
      pred.fLARS.test <- predict(res.fLARS, newldat)
    }
    MSEP.fLARS <- mean((matrix(pred.fLARS.test, ncol = 1) - Test_y)^2, na.rm = T)
    
    #############################################################################
    ## prediction using features selected by the reference paper SHFIM
    #############################################################################
    feature.SHFIM <- c(1, 2, 3, 5, 6, 7, 8)
    features.selected <- paste0("x", feature.SHFIM)
    res.SHFIM <- try(fregre.glm.vs(data=ldat, y="y", numbasis.opt=F
                                    , dcor.min=0, alpha = 1
                                    , include = features.selected
                                    # , trace = T
    ))
    if('try-error' %in% class(res.SHFIM)) {
      feature.SHFIM <- NULL
    }else{
      feature.SHFIM <- as.numeric(stringr::str_remove(res.SHFIM$ipredictors, pattern = "x"))
    }
    pred.SHFIM.train <- predict(res.SHFIM, ldat)
    feature.SHFIM
    # print(mean((matrix(pred.SHFIM.train, ncol = 1) - Train_y)^2, na.rm = T))
    
    # predict
    pred.SHFIM.test <- predict(res.SHFIM, newldat)
    MSEP.SHFIM <- mean((matrix(pred.SHFIM.test, ncol = 1) - Test_y)^2, na.rm = T)
    
    result.MSEP[irep, ] <- c(MSEP.fofcd.fb, MSEP.fofcd.osb, MSEP.fLARS, MSEP.glm, MSEP.SHFIM)
    print(result.MSEP[irep, ])
    result.feature[1:length(feature.fb), irep] <- feature.fb
    result.feature[8+(1:length(feature.osb)), irep] <- feature.osb
    result.feature[8*2+(1:length(feature.fLARS)), irep] <- feature.fLARS
    result.feature[8*3+(1:length(feature.glm)), irep] <- feature.glm
    result.feature[8*4+(1:length(feature.SHFIM)), irep] <- feature.SHFIM
  }
  colnames(result.MSEP) <- c("fofcd.fb", "fofcd.osb", "fLARS", "glm", "SHFIM")
  colnames(result.feature) <- paste0("rep_", 1:repitition)
  write.csv(result.MSEP, paste0("MSEP_GLM_of_", location, ".csv"))
  write.csv(result.feature, paste0("feature_of_", location, ".csv"))
}