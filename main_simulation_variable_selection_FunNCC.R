source("function_FOFCD.R")
library(fda.usc)
library(stringr)


iteration = 100
result.fb <- matrix(nrow = iteration, ncol = 20)
result.osb <- matrix(nrow = iteration, ncol = 20)

################################################################################
# variable selection 6.2 - functional complex nonlinear model
cat("variable selection 6.2 ", "functional complex nonlinear model\n")
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
  
  ## FunNCC
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  
  result.fb[iter, 1] <- FunNCC(y, list(X3_t), X = list(X1_t, X2_t), method = method, control = control)
  result.fb[iter, 2] <- FunNCC(y, list(X4_t), X = list(X1_t, X2_t), method = method, control = control)
  
  
  method = 'ospline_basis'; control = list(knots = c(tmin, tmin, seq(tmin, tmax, length.out = 5),
                                                     tmax, tmax),
                                           norder = 2, t = Time)
  result.osb[iter, 1] <- FunNCC(y, list(X3_t), X = list(X1_t, X2_t), method = method, control = control)
  result.osb[iter, 2] <- FunNCC(y, list(X4_t), X = list(X1_t, X2_t), method = method, control = control)
  
}

colnames(result.fb) <- colnames(result.osb) <- c("3|12", "4|12", 1:4, "2|1", "3|1", "4|1", "1|2", "3|2", "4|2",
                            "1|3", "2|3", "4|3", "1|4", "2|4", "3|4",
                            "12", "34")
write.csv(result.fb, paste0("result_variable_selection_simu_iteration2_FunNCC_fb.csv"), row.names = F)
write.csv(result.osb, paste0("result_variable_selection_selection_simu_iteration2_FunNCC_osb.csv"), row.names = F)
