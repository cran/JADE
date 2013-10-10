k_JADE <- function(X, k=1, eps = 1e-06, maxiter = 100, na.action = na.fail)
{
  X <- na.action(X)
  #if (!all(sapply(X, is.numeric)))  stop("'X' must be numeric")
  
  X <- as.matrix(X)
 
  n <- nrow(X)
  p <- ncol(X)   

  X<-scale(X,scale=F)
  Col.center <- attr(X,"scaled:center")

  W0 <- FOBI(X)$W
 
  Z <- X%*%t(W0) 

  CM <- array(0,c(p,p,sum((p-k+1):p)))
  R <- diag(p)
  Qij <- matrix(0,p,p)
  Zim <- numeric(p)
  Zjm <- numeric(p)
  scale2 <- rep(1,p)/n
  l<-1

  for(im in 1:p){
   Zim<-Z[,im]
   Qij<-t((Zim *Zim) %*% t(scale2) * Z) %*% Z - R-2*R[,im]%*%t(R[,im])
   CM[,,l]=Qij
   l<-l+1
   if(k>1){
    for(jm in 1:(k-1)){ 
     if(im+jm<=p){
      Zjm<-Z[,im+jm]
      Qij<-t((Zim *Zjm) %*% t(scale2) * Z) %*% Z - R[,im]%*%t(R[,im+jm])- 
            R[,im+jm]%*%t(R[,im])
      CM[,,l]=sqrt(2)*Qij
      l<-l+1   
     }
    }
   }
  }
  V<-rjd.fortran(CM,maxiter=maxiter,eps=eps)$V
 
  W <- t(V)%*%W0
  S <- X %*% t(W)
  A <- solve(W)
  colnames(S) <- paste("IC.", 1:p, sep="")
  res<-list(A=A, W=W, S=S, Xmu=Col.center)
  class(res) <- "bss"
  return(res)
}

