library(flars)
source("function_FOFCD.R")
source("function_FCCA-DC.R")


iteration = 100

################################################################################
# conditional correlation 1

n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)

col_simu <- c('FunNCC12.fb', 'FunNCC21.fb', 'FunNCC1.fb', 'FunNCC2.fb'
              ,'FunNCC12.raw', 'FunNCC21.raw', 'FunNCC1.raw', 'FunNCC2.raw' 
              ,'FunNCC12.osb', 'FunNCC21.osb', 'FunNCC1.osb', 'FunNCC2.osb'
              ,'FCCA12', 'FCCA1', 'FCCA2'
              , 'FCCA.1', 'FCCA.2','FSCA.1', 'FSCA.2'
              )
result.corr <- data.frame(matrix(ncol = length(col_simu), nrow = 0))
result.time <- data.frame(matrix(ncol = length(col_simu), nrow = 0))

for (iter in 1:iteration){
  cat("iteration times ", iter, "\n")
  x1 <- matrix(0, nrow = n, ncol = nt)
  x2 <- matrix(0, nrow = n, ncol = nt)
  w <- runif(n, min = 0, max = 1)
  v <- runif(n, min = 0, max = 1)
  for (idex in 1:n){
    x1[idex, ] <- sin(10*pi*w[idex]*Time)
    x2[idex, ] <- sin(10*pi*w[idex]*Time) + sin(5*pi*v[idex]*Time)
  }
  y <- rowSums(x2)/nt
  
  ## FunNCC
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  t1=Sys.time()
  FunNCC21.fb.corr <- FunNCC(y, list(x2=x2), list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC21.fb.time <- t2-t1
  t1=Sys.time()
  FunNCC12.fb.corr <- FunNCC(y, list(x1=x1), list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC12.fb.time <- t2-t1
  t1=Sys.time()
  FunNCC2.fb.corr <- FunNCC(y, list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC2.fb.time <- t2-t1
  t1=Sys.time()
  FunNCC1.fb.corr <- FunNCC(y, list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC1.fb.time <- t2-t1
  
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  t1=Sys.time()
  FunNCC21.osb.corr <- FunNCC(y, list(x2=x2), list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC21.osb.time <- t2-t1
  t1=Sys.time()
  FunNCC12.osb.corr <- FunNCC(y, list(x1=x1), list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC12.osb.time <- t2-t1
  t1=Sys.time()
  FunNCC2.osb.corr <- FunNCC(y, list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC2.osb.time <- t2-t1
  t1=Sys.time()
  FunNCC1.osb.corr <- FunNCC(y, list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC1.osb.time <- t2-t1
  
  # FCCA
  t1=Sys.time()
  FCCA1.corr <- flars::fccaGen(list(x1=x1), c(y), type='cor', method='basis')
  t2=Sys.time()
  FCCA1.time <- t2-t1
  t1=Sys.time()
  FCCA2.corr <- flars::fccaGen(list(x2=x2), c(y), type='cor', method='basis') 
  t2=Sys.time()
  FCCA2.time <- t2-t1
  t1=Sys.time()
  FCCA12.corr <- flars::fccaGen(list(x1=x1, x2=x2), y, type='cor', method='basis')
  t2=Sys.time()
  FCCA12.time <- t2-t1
  
  ## FSCA and FSCA-DC
  t1=Sys.time()
  FCCA.1.corr <- fcca.vecy(X=t(x1), Y=y, grid=Time, J = 5)$dcor.cca
  t2=Sys.time()
  FCCA.1.time <- t2-t1
  t1=Sys.time()
  FSCA1.corr <- fsca.vecy(X=t(x1), Y=y, grid=Time, J = 5)$dcor.sca
  t2=Sys.time()
  FSCA1.time <- t2-t1
  t1=Sys.time()
  FCCA.2.corr <- fcca.vecy(X=t(x2), Y=y, grid=Time, J = 5)$dcor.cca
  t2=Sys.time()
  FCCA.2.time <- t2-t1
  t1=Sys.time()
  FSCA2.corr <- fsca.vecy(X=t(x2), Y=y, grid=Time, J = 5)$dcor.sca
  t2=Sys.time()
  FSCA2.time <- t2-t1
  
  result.corr <- rbind(result.corr, c(FunNCC12.fb.corr, FunNCC21.fb.corr, FunNCC1.fb.corr, FunNCC2.fb.corr  
                                      ,FunNCC12.osb.corr, FunNCC21.osb.corr, FunNCC1.osb.corr, FunNCC2.osb.corr 
                                      ,FCCA12.corr, FCCA1.corr, FCCA2.corr
                                      , FCCA.1.corr, FCCA.2.corr, FSCA1.corr, FSCA2.corr
  ))
  result.time <- rbind(result.time, c(FunNCC12.fb.time, FunNCC21.fb.time, FunNCC1.fb.time, FunNCC2.fb.time 
                                      ,FunNCC12.osb.time, FunNCC21.osb.time, FunNCC1.osb.time, FunNCC2.osb.time 
                                      ,FCCA12.time, FCCA1.time, FCCA2.time 
                                      , FCCA.1.time, FCCA.2.time, FSCA1.time, FSCA2.time
  ))
}
colnames(result.corr) <- col_simu
colnames(result.time) <- col_simu
write.csv(result.corr, "result_conditional_correlation_simu_redundant_corr.csv", row.names = F)
write.csv(result.time, "result_conditional_correlation_simu_redundant_time.csv", row.names = F)



################################################################################
# conditional correlation 2
n=1000
Til=200
TRange=c(0, 1)
tmax = max(TRange)
tmin = min(TRange)
Time <- seq(from = tmin, to = tmax, length.out = Til)
nt <- length(Time)

col_simu <- c('FunNCC12.fb', 'FunNCC21.fb', 'FunNCC1.fb', 'FunNCC2.fb'
              ,'FunNCC12.raw', 'FunNCC21.raw', 'FunNCC1.raw', 'FunNCC2.raw' 
              ,'FunNCC12.osb', 'FunNCC21.osb', 'FunNCC1.osb', 'FunNCC2.osb'
              ,'FCCA12', 'FCCA1', 'FCCA2'
              , 'FCCA.1', 'FCCA.2','FSCA.1', 'FSCA.2'
              )
result.corr <- data.frame(matrix(ncol = length(col_simu), nrow = 0))
result.time <- data.frame(matrix(ncol = length(col_simu), nrow = 0))

for (iter in 1:iteration){
  cat("iteration times ", iter, "\n")
  w <- runif(n, min = 0, max = 1)
  v <- runif(n, min = 0, max = 1)
  eps <- rnorm(n, mean = 0, sd = 1)
  x1 <- matrix(0, nrow = n, ncol = nt)
  x2 <- matrix(0, nrow = n, ncol = nt)
  y <- array(0, dim = n)
  for (idex in 1:n){
    x1[idex, ] <- sin(10*pi*w[idex]*Time)
    x2[idex, ] <- sin(5*pi*v[idex]*Time)
  }
  y <- (w + v) %% 1
  
  ## FunNCC
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  t1=Sys.time()
  FunNCC21.fb.corr <- FunNCC(y, list(x2=x2), list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC21.fb.time <- t2-t1
  t1=Sys.time()
  FunNCC12.fb.corr <- FunNCC(y, list(x1=x1), list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC12.fb.time <- t2-t1
  t1=Sys.time()
  FunNCC2.fb.corr <- FunNCC(y, list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC2.fb.time <- t2-t1
  t1=Sys.time()
  FunNCC1.fb.corr <- FunNCC(y, list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC1.fb.time <- t2-t1
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  t1=Sys.time()
  FunNCC21.osb.corr <- FunNCC(y, list(x2=x2), list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC21.osb.time <- t2-t1
  t1=Sys.time()
  FunNCC12.osb.corr <- FunNCC(y, list(x1=x1), list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC12.osb.time <- t2-t1
  t1=Sys.time()
  FunNCC2.osb.corr <- FunNCC(y, list(x2=x2), method=method, control=control)
  t2=Sys.time()
  FunNCC2.osb.time <- t2-t1
  t1=Sys.time()
  FunNCC1.osb.corr <- FunNCC(y, list(x1=x1), method=method, control=control)
  t2=Sys.time()
  FunNCC1.osb.time <- t2-t1
  
  ## FCCA
  t1=Sys.time()
  FCCA1.corr <- flars::fccaGen(list(x1=x1), c(y), type='cor', method='basis')
  t2=Sys.time()
  FCCA1.time <- t2-t1
  t1=Sys.time()
  FCCA2.corr <- flars::fccaGen(list(x2=x2), c(y), type='cor', method='basis') 
  t2=Sys.time()
  FCCA2.time <- t2-t1
  t1=Sys.time()
  FCCA12.corr <- flars::fccaGen(list(x1=x1, x2=x2), y, type='cor', method='basis')
  t2=Sys.time()
  FCCA12.time <- t2-t1
  
  ## FSCA and FSCA-DC
  t1=Sys.time()
  FCCA.1.corr <- fcca.vecy(X=t(x1), Y=y, grid=Time, J = 5)$dcor.cca
  t2=Sys.time()
  FCCA.1.time <- t2-t1
  t1=Sys.time()
  FSCA1.corr <- fsca.vecy(X=t(x1), Y=y, grid=Time, J = 5)$dcor.sca
  t2=Sys.time()
  FSCA1.time <- t2-t1
  t1=Sys.time()
  FCCA.2.corr <- fcca.vecy(X=t(x2), Y=y, grid=Time, J = 5)$dcor.cca
  t2=Sys.time()
  FCCA.2.time <- t2-t1
  t1=Sys.time()
  FSCA2.corr <- fsca.vecy(X=t(x2), Y=y, grid=Time, J = 5)$dcor.sca
  t2=Sys.time()
  FSCA2.time <- t2-t1
  
  result.corr <- rbind(result.corr, c(FunNCC12.fb.corr, FunNCC21.fb.corr, FunNCC1.fb.corr, FunNCC2.fb.corr 
                                      ,FunNCC12.osb.corr, FunNCC21.osb.corr, FunNCC1.osb.corr, FunNCC2.osb.corr 
                                      ,FCCA12.corr, FCCA1.corr, FCCA2.corr
                                      , FCCA.1.corr, FCCA.2.corr, FSCA1.corr, FSCA2.corr
                                      ))
  result.time <- rbind(result.time, c(FunNCC12.fb.time, FunNCC21.fb.time, FunNCC1.fb.time, FunNCC2.fb.time 
                                      ,FunNCC12.osb.time, FunNCC21.osb.time, FunNCC1.osb.time, FunNCC2.osb.time 
                                      ,FCCA12.time, FCCA1.time, FCCA2.time 
                                      , FCCA.1.time, FCCA.2.time, FSCA1.time, FSCA2.time
                                      ))
  }
colnames(result.corr) <- col_simu
colnames(result.time) <- col_simu
write.csv(result.corr, "result_conditional_correlation_simu_corelevant_corr.csv", row.names = F)
write.csv(result.time, "result_conditional_correlation_simu_corelevant_time.csv", row.names = F)

