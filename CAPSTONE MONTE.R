library(MASS)
library(np)
###FUNCTIONS

kereg=function(y, Xdata, newX, h){                # y-data is in y and x-data in Xdata. newX is a single d-dim x-point to be classified.
  nn=length(y)
  dum1=matrix(rep(newX,nn), nrow=nn, byrow=TRUE)
  dum2=(Xdata-dum1)^2
  dum3=apply(dum2, 1, sum)/h^2
  answer=sum(y*exp(-dum3))
  return(answer) }

kern.reg=function(y, Xdata, X.points, h){         # X.points has more than one X.points
  n0=length(X.points[,1])
  answer=rep(0,n0)
  for(i in 1:n0){
    answer[i]= kereg(y, Xdata, X.points[i,], h)       # makes a call on the above function
  }
  return(answer) }

kern.reg_CROSS=function(y, Xdata, X.points, h){         # X.points has more than one X.points
  n0=length(X.points)
  answer=rep(0,n0)
  for(i in 1:n0){
    answer[i]= kereg(y, Xdata, X.points, h)       # makes a call on the above function
  }
  return(answer) }

CovMatrix=function(d) {           ##Covariance Matrix of Normal Vector
  SIGMA=matrix(0,nrow=d,ncol=d)
  for(i in 1:d){
    for(j in 1:d){
      SIGMA[i,j]=2^(-abs(i-j))
    }
  }
  return(SIGMA)
}

Create_Cauchy=function(n,d) {
  x=mvrnorm((n/2), mu=rep(1,d), Sigma=CovMatrix(d))
  x2=matrix(rcauchy((n/2)*d), (n/2), d)
  distX<-rbind(x,x2)
  return(distX)         
}

Create_Norm1=function(n,d) {
  x=mvrnorm((n/2), mu=rep(1,d), Sigma=CovMatrix(d))
  x2=matrix(runif(((n/2)*d), min=0, max=2), (n/2), d)
  distX<-rbind(x,x2)
  return(distX)
}

Create_Norm2=function(n,d) {
  x=mvrnorm((n/2), mu=rep(2,d), Sigma=CovMatrix(d))
  x2=matrix(runif(((n/2)*d), min=0, max=2), (n/2), d)
  distX<-rbind(x,x2)
  return(distX)
}
Create_Test_Cauchy=function(d, h_in) {
      newdistX <- Create_Cauchy(10000,d)
      newdistY <- c(rep(1,5000), rep(0, 5000))
      mykern <- kern.reg(distY, distX, newdistX, h=h_in)
      mykern2 <- kern.reg((1-distY), distX, newdistX, h=h_in)
      junk <- 1 * (mykern > mykern2)
      junk2 <- 1 * (junk == newdistY)
      error <- 1 - (sum(junk2)/length(junk2))
      return(error)
}
Create_Test_Norm1=function(d, h_in) {
  newdistX <- Create_Norm1(10000,d)
  newdistY <- c(rep(1,5000), rep(0, 5000))
  mykern <- kern.reg(distY, distX, newdistX, h=h_in)
  mykern2 <- kern.reg((1-distY), distX, newdistX, h=h_in)
  junk <- 1 * (mykern > mykern2)
  junk2 <- 1 * (junk == newdistY)
  error <- 1 - (sum(junk2)/length(junk2))
  return(error)
}
Create_Test_Norm2=function(d, h_in) {
  newdistX <- Create_Norm2(10000,d)
  newdistY <- c(rep(1,5000), rep(0, 5000))
  mykern <- kern.reg(distY, distX, newdistX, h=h_in)
  mykern2 <- kern.reg((1-distY), distX, newdistX, h=h_in)
  junk <- 1 * (mykern > mykern2)
  junk2 <- 1 * (junk == newdistY)
  error <- 1 - (sum(junk2)/length(junk2))
  return(error)
}

Find_holdout=function(distX, distY) {
 
  avg_error <- matrix(,nrow=200,ncol=20)
  
  for(j in 1:20) {
    ind = sample(1:nrow(distX), nrow(distX) * .7)
    trainX <- distX[ind,]
    testX <- distX[-ind,]
    trainY <- distY[ind]
    testY <- distY[-ind]
    
    h <- .0001
    h_vec <- NULL
    error_vec <- NULL
    
    for(i in 1:200) {
      mykern <- kern.reg(trainY, trainX, testX, h)
      mykern2 <- kern.reg((1-trainY), trainX, testX, h)
      junk <- 1 * (mykern > mykern2)
      junk2 <- 1 * (junk == testY)
      error <- 1 - (sum(junk2)/length(junk2))
      h_vec <- c(h_vec, h)
      error_vec <- c(error_vec, error)
      h = h + .01
    }
    avg_error[,j] = error_vec
  }
  
  hold_avg <- rowSums(avg_error)/20
  hold_h=h_vec[which.min(hold_avg)]
  return(hold_h)
}

Find_cross=function(distX, distY,d) {
  h <- .0001
  h_vec <- NULL
  error_vec <- NULL
  mykern <- NULL
  mykern2 <- NULL
  hold_y <- NULL
  
  for(j in 1:200) {
    for (i in 1:length(distY)){
      holdX = distX[i,]
      holdY = distY[i]
      this_dataX = distX[-i,]
      this_dataY = distY[-i]
      mykern <- c(mykern, kern.reg_CROSS(this_dataY, this_dataX, holdX, h))
      mykern2 <- c(mykern2, kern.reg_CROSS((1 - this_dataY), this_dataX, holdX, h))
      hold_y <- c(hold_y, rep(holdY,d))
      
    }
    junk <- 1 * (mykern > mykern2)
    junk2 <- 1 * (junk == hold_y)
    error <- 1 - (sum(junk2)/length(junk2))
    h_vec <- c(h_vec, h)
    error_vec <- c(error_vec, error)
    h = h + .01
    mykern <- NULL
    mykern2 <- NULL
    hold_y <- NULL
  }
  
  cross_h=h_vec[which.min(error_vec)]
  return(cross_h)
}

Find_npreg = function(distX, distY,d) {
  deleted <- npregbw(xdat = distX, ydat = distY)
  vals <- 1 * deleted$bandwidth$x
  npreg_h <- prod(vals)^(1/d)
  return(npreg_h)
}

h_vec_cross <- NULL
h_vec_hold <- NULL
h_vec_npreg <- NULL
err_cross <- NULL
err_hold <- NULL
err_npreg <- NULL
err_vec_cross <- NULL
err_vec_hold <- NULL
err_vec_npreg <- NULL

n=50
d=10 

for (i in 1:50) {
  distX <- Create_Cauchy(n,d)
  distY <- c(rep(1,(n/2)),rep(0,(n/2)))
  cross_h <- Find_cross(distX, distY,d)
  hold_h <- Find_holdout(distX, distY)
  npreg_h <- Find_npreg(distX, distY,d)
  h_vec_cross <- c(h_vec_cross, cross_h)
  h_vec_hold <- c(h_vec_hold, hold_h)
  h_vec_npreg <- c(h_vec_npreg, npreg_h)
  err_cross <- Create_Test_Cauchy(d,cross_h)
  err_hold <- Create_Test_Cauchy(d,hold_h)
  err_npreg <- Create_Test_Cauchy(d,npreg_h)
  err_vec_cross <- c(err_vec_cross, err_cross)
  err_vec_hold <- c(err_vec_hold, err_hold)
  err_vec_npreg <- c(err_vec_npreg, err_npreg)
}
for (i in 1:238) {
  err_cross <- Create_Test_Cauchy(10,h_vec_cross[i])
  err_hold <- Create_Test_Cauchy(10,h_vec_hold[i])
  err_npreg <- Create_Test_Cauchy(10,h_vec_npreg[i])
  err_vec_cross <- c(err_vec_cross, err_cross)
  err_vec_hold <- c(err_vec_hold, err_hold)
  err_vec_npreg <- c(err_vec_npreg, err_npreg)
}

distX <- Create_Norm1(20,10)
distY <- c(rep(1,10),rep(0,10))
cross_h <- Find_cross(distX, distY,10)
hold_h <- Find_holdout(distX, distY)
npreg_h <- Find_npreg(distX, distY,10)
Create_Test_Norm1(10,cross_h)
Create_Test_Norm1(10,hold_h)
Create_Test_Norm1(10,npreg_h)

distX <- Create_Norm2(20,10)
distY <- c(rep(1,10),rep(0,10))
cross_h <- Find_cross(distX, distY,10)
hold_h <- Find_holdout(distX, distY)
npreg_h <- Find_npreg(distX, distY,10)
Create_Test_Norm2(10,cross_h)
Create_Test_Norm2(10,hold_h)
Create_Test_Norm2(10,npreg_h)

Simul_Cauchy = function(n,d) {
  h_vec_cross <- NULL
  h_vec_hold <- NULL
  h_vec_npreg <- NULL
  err_cross <- NULL
  err_hold <- NULL
  err_npreg <- NULL
  distX <- Create_Cauchy(n,d)
  distY <- c(rep(1,(n/2)),rep(0,(n/2)))
  cross_h <- Find_cross(distX, distY,d)
  hold_h <- Find_holdout(distX, distY)
  npreg_h <- Find_npreg(distX, distY,d)
  h_vec_cross <- c(h_vec_cross, cross_h)
  h_vec_hold <- c(h_vec_hold, hold_h)
  h_vec_npreg <- c(h_vec_npreg, npreg_h)
  err_cross <- Create_Test_Cauchy(d,cross_h)
  err_hold <- Create_Test_Cauchy(d,hold_h)
  err_npreg <- Create_Test_Cauchy(d,npreg_h)
  results <- matrix(nrow = 3, ncol = 2)
}

Simul_Cauchy(20,10)

n = 20
d = 10
distx <- Create_Cauchy(n,d)
distY<- c(rep(1,(n/2)),rep(0,(n/2)))
cross_h <- Find_cross(distX, distY,d)
h_vec_cross <- c(h_vec_cross, cross_h)
