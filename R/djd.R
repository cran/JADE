# Function for deflation based JD
#
# input: 
#   X = p,p,k array with k pxp matrices which should be diagonalized
#   G = for optimization, options are "pow" and "log", default is "pow"
#   r = power used, if G="pow", default is 2 
#   eps = convergence criterion
#   maxiter = maximum number of iterations
#
# output:
#   orthogonal matrix

djd <- function(X, G="pow", r=2, eps = 1e-06, maxiter = 100)
    {
    G <- match.arg(G, c("pow","log"))
    
    WN <- eigen(X[,,1])$vectors   
    p <- dim(X)[1]
    k <- dim(X)[3]
    W <- matrix(0,p,p)
    
    W <- switch(G,
        "pow"={
               djd.pow(X=X, W=W,r=r, WN=WN, k=k, p=p, eps=eps, maxiter=maxiter)
               }
        ,
        "log"={
              djd.log(X=X, W=W, WN=WN, k=k, p=p, eps=eps, maxiter=maxiter)
              }
        )
    
     W
     } 

########################################################
#
# djd subfunctions
#
djd.pow <- function(X,W,r,WN,k,p,eps,maxiter)
    {
    for (i in 1:p){
                    wn <- (WN[,i, drop=FALSE])
                    for (it in 1:maxiter){
                        w<-wn
                        wn<- matrix(0,p,1)
                        for (mi in 1:k){
                            wn <- wn +  r*(as.numeric(t(w) %*% X[,,mi] %*% w))^(r-1) * X[,,mi] %*% w}
                        wn <- wn- crossprod(W)%*%wn
                        wn <- wn / sqrt(sum(wn^2)) 
                        if (sqrt(sum((w-wn)^2))<eps || sqrt(sum((w+wn)^2))<eps) break
                        if (it==maxiter) stop("no convergence reached")
                        }
                    W[i,]<-t(wn)
                    }
    t(W)
    }
    
        
djd.log <- function(X,W,WN,k,p,eps,maxiter)
    {
    for (i in 1:p){
                    wn <- (WN[,i, drop=FALSE])
                    for (it in 1:maxiter){
                        w<-wn
                        wn<- matrix(0,p,1)
                        for (mi in 1:k){
                            wn <- wn +  1/as.numeric(crossprod(w, X[,,mi]) %*% w)* crossprod(X[,,mi], w)}
                        wn <- wn- crossprod(W)%*%wn
                        wn <- wn / sqrt(sum(wn^2)) 
                        if (sqrt(sum((w-wn)^2))<eps || sqrt(sum((w+wn)^2))<eps) break
                        if (it==maxiter) stop("no convergence reached")
                        }
                    W[i,]<-t(wn)
                    }
    t(W)
    }
