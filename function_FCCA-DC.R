require(fda)
require(magic)
require(mvtnorm)
require(splines)

CCA<-function(X,Y){
  Xsvd<-svd(t(X)%*%X); Ysvd<-svd(t(Y)%*%Y)
  
  if (dim(t(X))[1] == 1){
    Xinv<-Xsvd$u%*%(1/(Xsvd$d+1e-10))%*%t(Xsvd$u)
    Xinvhalf<-Xsvd$u%*%(1/sqrt(Xsvd$d+1e-10))%*%t(Xsvd$u)
  }else{
    Xinv<-Xsvd$u%*%diag(1/(Xsvd$d+1e-10))%*%t(Xsvd$u)
    Xinvhalf<-Xsvd$u%*%diag(1/sqrt(Xsvd$d+1e-10))%*%t(Xsvd$u)
  }    
  
  if (dim(t(Y))[1] == 1){
    Yinv<-Ysvd$u%*%(1/(Ysvd$d+1e-10))%*%t(Ysvd$u)
    Yinvhalf<-Ysvd$u%*%(1/sqrt(Ysvd$d+1e-10))%*%t(Ysvd$u)
  }else{
    Yinv<-Ysvd$u%*%diag(1/(Ysvd$d+1e-10))%*%t(Ysvd$u)
    Yinvhalf<-Ysvd$u%*%diag(1/sqrt(Ysvd$d+1e-10))%*%t(Ysvd$u)
  }
  temp<-Xinvhalf%*%t(X)%*%Y%*%Yinv%*%t(Y)%*%X%*%Xinvhalf
  alpha<-as.vector(Xinvhalf%*%svd(temp)$u[,1])
  
  temp<-Yinvhalf%*%t(Y)%*%X%*%Xinv%*%t(X)%*%Y%*%Yinvhalf
  beta<-as.vector(Yinvhalf%*%svd(temp)$u[,1])
  
  return(list(alpha=alpha/sqrt(sum(alpha^2)),beta=beta/sqrt(sum(beta^2))))
}

SCA<-function(X,Y){
  XYsvd<-svd(t(X)%*%Y)
  alpha<-XYsvd$u[,1]; beta<-XYsvd$v[,1]
  return(list(alpha=alpha/sqrt(sum(alpha^2)),beta=beta/sqrt(sum(beta^2))))
}

dCovCCA<-function(alpha,beta,X,Y){
  tstep<-0.001
  alpha<-alpha/sqrt(sum(alpha^2)); beta<-beta/sqrt(sum(beta^2))
  for (iter in 1:100){
    temp<-computedCovgrad(alpha,beta,X,Y)
    alpha<-alpha+tstep*temp$gradalpha
    beta<-beta+tstep*temp$gradbeta
    
    #normalize
    alpha<-alpha/sqrt(sum(alpha^2)); beta<-beta/sqrt(sum(beta^2));
    
    #cat(computedCor(alpha,beta,X,Y)[1:2],"\n")
  }
  return(list(alpha=alpha,beta=beta))
}




dCorCCA2<-function(alpha,beta,X,Y){
  tstep<-0.003
  for (iter in 1:100){
    
    temp<-computedCorgrad(alpha,beta,X,Y)
    alpha<-alpha+tstep*temp$gradalpha
    beta<-beta+tstep*temp$gradbeta
    
    alpha<-alpha/sqrt(sum(alpha^2)); beta<-beta/sqrt(sum(beta^2))
    
    #cat(computedCor(alpha,beta,X,Y)[1:2],"\n")
  }
  
  return(list(alpha=alpha,beta=beta))
}


computedCov<-function(alpha,beta,X,Y){
  if (length(alpha)==1) X <- as.matrix(X, ncol = 1)
  if (length(beta)==1) Y <- as.matrix(Y, ncol = 1)
  
  x<-as.vector(X%*%alpha)
  y<-as.vector(Y%*%beta)
  n<-length(x)
  tildeA<-abs(outer(x,x,"-"))
  tildeB<-abs(outer(y,y,"-"))
  H<-diag(n)-matrix(1,nrow=n,ncol=n)/n
  A<-H%*%tildeA%*%H; B<-H%*%tildeB%*%H;
  sum(A*B)/n^2
}

computedCor<-function(alpha,beta,X,Y){
  if (length(alpha)==1) X <- as.matrix(X, ncol = 1)
  if (length(beta)==1) Y <- as.matrix(Y, ncol = 1)
  
  x<-as.vector(X%*%alpha)
  y<-as.vector(Y%*%beta)
  n<-length(x)
  tildeA<-abs(outer(x,x,"-"))
  tildeB<-abs(outer(y,y,"-"))
  H<-diag(n)-matrix(1,nrow=n,ncol=n)/n
  A<-H%*%tildeA%*%H; B<-H%*%tildeB%*%H;
  return(c( sum(A*B)/sqrt(sum(A^2)*sum(B^2)), sum(A*B)/n/n, sqrt(sum(A^2))/n, sqrt(sum(B^2))/n ) )
}

computedCovgrad<-function(alpha,beta,X,Y){
  if (length(alpha)==1) X <- as.matrix(X, ncol = 1)
  if (length(beta)==1) Y <- as.matrix(Y, ncol = 1)
  
  x<-as.vector(X%*%alpha)
  y<-as.vector(Y%*%beta)
  n<-length(x)
  
  Cx<-(outer(x,x,"-")>0)*2-1
  Cy<-(outer(y,y,"-")>0)*2-1
  tildeA<-abs(outer(x,x,"-"))
  tildeB<-abs(outer(y,y,"-"))
  H<-diag(n)-matrix(1,nrow=n,ncol=n)/n
  A<-H%*%tildeA%*%H; B<-H%*%tildeB%*%H;
  
  D<-B*Cx; diag(D)<-0
  gradalpha<-t(apply(D,1,sum))%*%X
  D<-A*Cy; diag(D)<-0
  gradbeta<-t(apply(D,1,sum))%*%Y
  
  return(list(gradalpha=as.vector(gradalpha),gradbeta=as.vector(gradbeta)))
}


computedCorgrad<-function(alpha,beta,X,Y){
  if (length(alpha)==1) X <- as.matrix(X, ncol = 1)
  if (length(beta)==1) Y <- as.matrix(Y, ncol = 1)
  
  x<-as.vector(X%*%alpha)
  y<-as.vector(Y%*%beta)
  n<-length(x)
  
  Cx<-(outer(x,x,"-")>0)*2-1
  Cy<-(outer(y,y,"-")>0)*2-1
  tildeA<-abs(outer(x,x,"-"))
  tildeB<-abs(outer(y,y,"-"))
  H<-diag(n)-matrix(1,nrow=n,ncol=n)/n
  A<-H%*%tildeA%*%H; B<-H%*%tildeB%*%H;
  
  Dx<-B*Cx; diag(Dx)<-0
  Dy<-A*Cy; diag(Dy)<-0
  
  Da<-A*Cx; diag(Da)<-0
  Db<-B*Cy; diag(Db)<-0
  
  gradalpha<-t(apply(Dx,1,sum))%*%X/sum(A*B)-t(apply(Da,1,sum))%*%X/sum(A^2)
  
  gradbeta<-t(apply(Dy,1,sum))%*%Y/sum(A*B)-t(apply(Db,1,sum))%*%Y/sum(B^2)
  
  return(list(gradalpha=as.vector(gradalpha),gradbeta=as.vector(gradbeta)))
}


initialize.old<-function(X,Y,nsample=10){
  n<-dim(X)[1]
  curmax<-0
  curmaxi<-1; curmaxj<-1
  seqi<-sort(sample(1:n,nsample))
  seqj<-sort(sample(1:n,nsample))
  for (i in seqi){
    for (j in seqj){
      temp<-computedCor(X[i,],Y[j,],X,Y)[1]
      #cat(i,j,temp,"\n")
      if (temp>curmax){ curmax<-temp; curmaxi<-i; curmaxj<-j}
    }
  }
  return(list(alpha=as.vector(X[curmaxi,]),beta=as.vector(Y[curmaxj,])))
}



initialize<-function(X,Y,nsample=10){
  n<-dim(X)[1]
  curmax<-0
  seqi<-sort(sample(1:n,nsample))
  for (i in seqi){
    seqj<-sort(sample(1:n,nsample))
    for (j in seqj){
      temp<-computedCor(X[i,],Y[j,],X,Y)[1]
      if (temp>curmax){ curmax<-temp; curalpha<-X[i,]; curbeta<-Y[j,]}
    }
  }
  return(list(alpha=as.vector(curalpha),beta=as.vector(curbeta)))
}


fcca <- function(X, Y, grid, J = c(5, 5)){
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  oriY<-oriY-apply(oriY,1,mean) #oriX,Y never change
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  myfdy<-Data2fd(argvals=grid,y=funY,mybasis) # number column of X is n
  mypcay<-pca.fd(myfdy,nharm=mybasis$nbasis-1)
  myscoresy<-mypcay$score #now number of ROW is n
  
  allX<-myscoresx; allY<-myscoresy
  allX<-t(t(allX)-apply(allX,2,mean))
  allY<-t(t(allY)-apply(allY,2,mean))
  
  ###################################
  ####### CCA #######################
  ###################################
  cat("CCA\n")
  
  J1 <- J[1]
  J2 <- J[2]
  temp<-CCA(as.matrix(allX[,1:J1]),as.matrix(allY[,1:J2]))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J1]%*%alpha)
  betafun<-as.vector(mybasismatrix%*%mypcay$harmonics$coefs[,1:J2]%*%beta)
  alphafun<-alphafun/sqrt(sum(alphafun^2)); betafun<-betafun/sqrt(sum(betafun^2))
  
  
  (cor.cca<-cor(t(funX)%*%alphafun,t(funY)%*%betafun))
  (dcor.cca<-computedCor(alphafun,betafun,t(funX),t(funY))[1])
  
  
  op <- par(mfrow=c(1,3))
  plot(alphafun,xlab="",ylab="",main="alpha")
  plot(betafun,xlab="",ylab="",main="beta")
  plot(alphafun%*%funX,betafun%*%funY,xlab="",ylab="",main="")
  par(op)
  
  return(list(cor.cca = cor.cca, dcor.cca = dcor.cca))
  
}

fsca <- function(X, Y, grid, J = c(5, 5)){
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  oriY<-oriY-apply(oriY,1,mean) #oriX,Y never change
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  myfdy<-Data2fd(argvals=grid,y=funY,mybasis) # number column of X is n
  mypcay<-pca.fd(myfdy,nharm=mybasis$nbasis-1)
  myscoresy<-mypcay$score #now number of ROW is n
  
  allX<-myscoresx; allY<-myscoresy
  allX<-t(t(allX)-apply(allX,2,mean))
  allY<-t(t(allY)-apply(allY,2,mean))
  
  ###################################
  ####### SCA #######################
  ###################################
  cat("SCA\n")
  
  
  
  J1 <- J[1]
  J2 <- J[2]
  temp<-SCA(as.matrix(allX[,1:J1]),as.matrix(allY[,1:J2]))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J1]%*%alpha)
  betafun<-as.vector(mybasismatrix%*%mypcay$harmonics$coefs[,1:J2]%*%beta)
  alphafun<-alphafun/sqrt(sum(alphafun^2)); betafun<-betafun/sqrt(sum(betafun^2))
  
  
  (cor.sca<-cor(t(funX)%*%alphafun,t(funY)%*%betafun))
  (dcor.sca<-computedCor(alphafun,betafun,t(funX),t(funY))[1])
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(betafun,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,betafun%*%funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.sca = cor.sca, dcor.sca = dcor.sca))
  
}


fcca.dc <- function(X, Y, grid, J = c(5, 5)){
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  oriY<-oriY-apply(oriY,1,mean) #oriX,Y never change
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  myfdy<-Data2fd(argvals=grid,y=funY,mybasis) # number column of X is n
  mypcay<-pca.fd(myfdy,nharm=mybasis$nbasis-1)
  myscoresy<-mypcay$score #now number of ROW is n
  
  allX<-myscoresx; allY<-myscoresy
  allX<-t(t(allX)-apply(allX,2,mean))
  allY<-t(t(allY)-apply(allY,2,mean))
  
  ##############################
  ####### dCor-CCA  ###########
  ##############################
  cat("dCor-CCA\n")
  
  J1 <- J[1]
  J2 <- J[2]
  init<-initialize(allX[,1:J1],allY[,1:J2])
  initalpha<-init$alpha;initbeta<-init$beta
  
  temp<-dCorCCA2(initalpha[1:J1],initbeta[1:J2],as.matrix(allX[,1:J1]),as.matrix(allY[,1:J2]))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J1]%*%alpha)
  betafun<-as.vector(mybasismatrix%*%mypcay$harmonics$coefs[,1:J2]%*%beta)
  alphafun<-alphafun/sqrt(sum(alphafun^2)); betafun<-betafun/sqrt(sum(betafun^2))
  
  (cor.dcor<-cor(t(funX)%*%alphafun,t(funY)%*%betafun))
  (dcor.dcor<-computedCor(alphafun,betafun,t(funX),t(funY))[1])
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(betafun,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,betafun%*%funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.dcor = cor.dcor, dcor.dcor = dcor.dcor))
}


fsca.dc <- function(X, Y, grid, J = c(5, 5)){
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  oriY<-oriY-apply(oriY,1,mean) #oriX,Y never change
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  myfdy<-Data2fd(argvals=grid,y=funY,mybasis) # number column of X is n
  mypcay<-pca.fd(myfdy,nharm=mybasis$nbasis-1)
  myscoresy<-mypcay$score #now number of ROW is n
  
  allX<-myscoresx; allY<-myscoresy
  allX<-t(t(allX)-apply(allX,2,mean))
  allY<-t(t(allY)-apply(allY,2,mean))
  
  ##############################
  ####### dCov-CCA #############
  ##############################
  cat("dCoV-CCA\n")
  
  J1 <- J[1]
  J2 <- J[2]
  init<-initialize(allX[,1:J1],allY[,1:J2])
  initalpha<-init$alpha;initbeta<-init$beta
  
  temp<-dCovCCA(initalpha[1:J1],initbeta[1:J2],as.matrix(allX[,1:J1]),as.matrix(allY[,1:J2]))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J1]%*%alpha)
  betafun<-as.vector(mybasismatrix%*%mypcay$harmonics$coefs[,1:J2]%*%beta)
  alphafun<-alphafun/sqrt(sum(alphafun^2)); betafun<-betafun/sqrt(sum(betafun^2))
  
  
  (cor.dcov<-cor(t(funX)%*%alphafun,t(funY)%*%betafun))
  (dcor.dcov<-computedCor(alphafun,betafun,t(funX),t(funY))[1])
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(betafun,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,betafun%*%funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.dcov = cor.dcov, dcor.dcov = dcor.dcov))
}


fcca.vecy <- function(X, Y, grid, J = 5){
  Y <- as.matrix(Y, ncol = 1)
  
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n

  allX<-myscoresx
  allX<-t(t(allX)-apply(allX,2,mean))
  
  ###################################
  ####### CCA #######################
  ###################################
  cat("CCA\n")
  
  temp<-CCA(as.matrix(allX[,1:J]),as.matrix(funY))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J]%*%alpha)
  alphafun<-alphafun/sqrt(sum(alphafun^2))
  
  
  (cor.cca<-cor(t(funX)%*%alphafun,beta*funY))
  (dcor.cca<-computedCor(alphafun,beta,t(funX),funY)[1])
  
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(beta,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,beta*funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.cca = cor.cca, dcor.cca = dcor.cca))
}

fsca.vecy <- function(X, Y, grid, J = 5){
  Y <- as.matrix(Y, ncol = 1)
  
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  allX<-myscoresx
  allX<-t(t(allX)-apply(allX,2,mean))
  
  ###################################
  ####### SCA #######################
  ###################################
  cat("SCA\n")
  
  temp<-SCA(as.matrix(allX[,1:J]),as.matrix(funY))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J]%*%alpha)
  alphafun<-alphafun/sqrt(sum(alphafun^2))
  
  
  (cor.sca<-cor(t(funX)%*%alphafun,beta*funY))
  (dcor.sca<-computedCor(alphafun,beta,t(funX),funY)[1])
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(beta,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,beta*funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.sca = cor.sca, dcor.sca = dcor.sca))
}


fcca.dc.vecy <- function(X, Y, grid, J = 5){
  Y <- as.matrix(Y, ncol = 1)
  
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  allX<-myscoresx
  allX<-t(t(allX)-apply(allX,2,mean))
  
  ##############################
  ####### dCor-CCA  ###########
  ##############################
  cat("dCor-CCA\n")
  
  init<-initialize(allX[,1:J],funY)
  initalpha<-init$alpha;initbeta<-init$beta
  
  temp<-dCorCCA2(initalpha[1:J],initbeta,as.matrix(allX[,1:J]),as.matrix(funY))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J]%*%alpha)
  alphafun<-alphafun/sqrt(sum(alphafun^2))
  
  (cor.dcor<-cor(t(funX)%*%alphafun,beta*funY))
  (dcor.dcor<-computedCor(alphafun,beta,t(funX),funY)[1])
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(beta,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,beta*funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.dcor = cor.dcor, dcor.dcor = dcor.dcor))
}


fsca.dc.vecy <- function(X, Y, grid, J = 5){
  Y <- as.matrix(Y, ncol = 1)
  
  ninternal=16 # number basis is ninternal+4
  mybreaks<-seq(range(grid)[1],range(grid)[2],length=ninternal+2) 
  mybasis<-create.bspline.basis(rangeval=range(grid),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(grid,mybasis) #ntimes*nbasis matrix
  
  oriX<-X; oriY<-Y
  oriX<-oriX-apply(oriX,1,mean)
  funX<-oriX; funY<-oriY 
  
  #use FPCA to project the functional data to functional scores
  myfdx<-Data2fd(argvals=grid,y=funX,mybasis) # number column of X is n
  mypcax<-pca.fd(myfdx,nharm=mybasis$nbasis-1)
  myscoresx<-mypcax$score #now number of ROW is n
  
  allX<-myscoresx
  allX<-t(t(allX)-apply(allX,2,mean))
  
  ##############################
  ####### dCov-CCA #############
  ##############################
  cat("dCoV-CCA\n")
  
  init<-initialize(allX[,1:J],funY)
  initalpha<-init$alpha;initbeta<-init$beta
  
  temp<-dCovCCA(initalpha[1:J],initbeta,as.matrix(allX[,1:J]),as.matrix(funY))
  alpha<-temp$alpha; beta<-temp$beta
  alphafun<-as.vector(mybasismatrix%*%mypcax$harmonics$coefs[,1:J]%*%alpha)
  alphafun<-alphafun/sqrt(sum(alphafun^2))
  
  
  (cor.dcov<-cor(t(funX)%*%alphafun,beta*funY))
  (dcor.dcov<-computedCor(alphafun,beta,t(funX),funY)[1])
  
  # op <- par(mfrow=c(1,3))
  # plot(alphafun,xlab="",ylab="",main="alpha")
  # plot(beta,xlab="",ylab="",main="beta")
  # plot(alphafun%*%funX,beta*funY,xlab="",ylab="",main="")
  # par(op)
  
  return(list(cor.dcov = cor.dcov, dcor.dcov = dcor.dcov))
}
