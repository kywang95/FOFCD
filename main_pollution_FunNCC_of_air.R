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
location <- "Guanyuan"
data <- read.csv(paste0(path, "PRSA_Data_", location, "_20130301-20170228.csv"), header = T)
data.2015 <- data[which(data$year=="2015"), ]

n <- nrow(data.2015)/24
# define the time points on which the functional predictor is observed. 
timepts <- unique(data.2015$hour)/24
ntime <- length(timepts)
Trange <- range(timepts)

repitition <- 10
result.fb <- matrix(nrow = repitition, ncol = 10)
result.osb <- matrix(nrow = repitition, ncol = 10)
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
  
  ntrain <- 300
  Train <- lapply(data, function(x) x[1:ntrain,])
  # Train <- data
  
  Train_x <- Train[-1]
  Train_y <- matrix(Train[[1]], ncol = 1)
  colnames(Train_y) <- names(data)[1]
  
  ndimX = lapply(Train_x, ncol)
  names(Train_x) <- names(data)[-1]
  
  # define the fourier basis 
  nbasis = 15
  temp_basis = create.fourier.basis(Trange, nbasis)
  # convert the functional predictor in training data into fda objects
  p <- length(Train_x)
  train_n <- nrow(Train_x[[1]])
  train_fda = array(dim = c(nbasis, train_n, p))
  train_data <- list()
  for(ip in 1:p){
    train_x_data =  Data2fd(timepts, t(Train_x[[ip]]), temp_basis)
    train_data[[ip]] <- t(eval.fd(timepts, train_x_data))
    func_cov = train_x_data$coefs
    train_fda[,,ip] <- func_cov
  }
  names(train_data) <- paste0("X", c(1:p))
  
  ## FunNCC
  method = 'fourier_basis'; control = list(nbasis = 2, norder = 2, t = timepts)
  result.fb[irep, 1] <- FunNCC(Train_y, train_data[4], X = train_data[c(1,6,7)], method = method, control = control)
  result.fb[irep, 2] <- FunNCC(Train_y, train_data[c(2,3)], X = train_data[c(1,6,7)], method = method, control = control)
  result.fb[irep, 3] <- FunNCC(Train_y, train_data[4], X = NULL, method = method, control = control)
  result.fb[irep, 4] <- FunNCC(Train_y, train_data[4], X = train_data[c(2,3)], method = method, control = control)
  result.fb[irep, 5] <- FunNCC(Train_y, train_data[2], X = NULL, method = method, control = control)
  result.fb[irep, 6] <- FunNCC(Train_y, train_data[2], X = train_data[4], method = method, control = control)
  result.fb[irep, 7] <- FunNCC(Train_y, train_data[3], X = NULL, method = method, control = control)
  result.fb[irep, 8] <- FunNCC(Train_y, train_data[3], X = train_data[4], method = method, control = control)
  result.fb[irep, 9] <- FunNCC(Train_y, train_data[c(2,3)], X = NULL, method = method, control = control)
  result.fb[irep, 10] <- FunNCC(Train_y, train_data[c(2,3)], X = train_data[4], method = method, control = control)
  
  
  tmax = max(timepts)
  tmin = min(timepts)
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = timepts)
  result.osb[irep, 1] <- FunNCC(Train_y, train_data[4], X = train_data[c(1,6,7)], method = method, control = control)
  result.osb[irep, 2] <- FunNCC(Train_y, train_data[c(2,3)], X = train_data[c(1,6,7)], method = method, control = control)
  result.osb[irep, 3] <- FunNCC(Train_y, train_data[4], X = NULL, method = method, control = control)
  result.osb[irep, 4] <- FunNCC(Train_y, train_data[4], X = train_data[c(2,3)], method = method, control = control)
  result.osb[irep, 5] <- FunNCC(Train_y, train_data[2], X = NULL, method = method, control = control)
  result.osb[irep, 6] <- FunNCC(Train_y, train_data[2], X = train_data[4], method = method, control = control)
  result.osb[irep, 7] <- FunNCC(Train_y, train_data[3], X = NULL, method = method, control = control)
  result.osb[irep, 8] <- FunNCC(Train_y, train_data[3], X = train_data[4], method = method, control = control)
  result.osb[irep, 9] <- FunNCC(Train_y, train_data[c(2,3)], X = NULL, method = method, control = control)
  result.osb[irep, 10] <- FunNCC(Train_y, train_data[c(2,3)], X = train_data[4], method = method, control = control)
}
colnames(result.fb) <- c("CO|167", "SO2NO2|167", "CO", "CO|SO2NO2", "SO2", "SO2|CO", "NO2", "NO2|CO", "SO2NO2", "SO2NO2|CO")
colnames(result.osb) <- c("CO|167", "SO2NO2|167", "CO", "CO|SO2NO2", "SO2", "SO2|CO", "NO2", "NO2|CO", "SO2NO2", "SO2NO2|CO")

write.csv(result.fb, paste0("FunNCC_fb_of_", location, ".csv"))
write.csv(result.osb, paste0("FunNCC_osb_of_", location, ".csv"))
