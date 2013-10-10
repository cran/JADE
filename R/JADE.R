`JADE` <-
function(X,n.comp=NULL, eps = 1e-06, maxiter = 100, na.action = na.fail)
    {
    # default eps would be 1/sqrt(N)/100
    X <- na.action(X)
    if (!all(sapply(X, is.numeric)))  stop("'X' must be numeric")
    X <- as.matrix(X)
    data.matrix<-X
    N <- dim(X)[1]
    X.cols <- dim(X)[2]
    
    if (is.null(n.comp)) n.comp <- X.cols
    if (n.comp > X.cols) stop("'n.comp' must be smaller than number of columns of X")
    
    X<-scale(X,scale=F)
    Col.center <- attr(X,"scaled:center")
    data.X<-X
    eigen.X<-eigen(t(X)%*%X/N)
    U1<-eigen.X$vectors
    D1<-eigen.X$values
    
    puiss<-sort(D1,decreasing=F)
    k<-1:length(D1)
    U<-U1[,rev(k)]
    rangeW <- (X.cols-n.comp+1):X.cols
    scales<-sqrt(puiss[rangeW])
    
    if (n.comp>1)
    {
    W <- diag(1/scales) %*% t(U[,rangeW])
    iW <- U[,rangeW] %*% diag(scales)
    }
    else
    {
    W<- t(U[,rangeW])/scales
    iW <- U[,rangeW]*scales
    }

    X <- X %*% t(W)

    dimsymm <- (n.comp*(n.comp+1))/2
    nbcm  <- dimsymm
    CM <- matrix(0,nrow=n.comp*nbcm,ncol=n.comp)
    R <- diag(n.comp)
    Qij <- matrix(0,ncol=n.comp,nrow=n.comp)
    Xim <- numeric(n.comp)
    Xjm <- numeric(n.comp)
    scale2 <- rep(1,n.comp)/N

    Range <- 1:n.comp

    for (im in 1:n.comp)
    {
    Xim<-X[,im]
    Qij<-t((Xim *Xim) %*% t(scale2) * X) %*% X - R-2*R[,im]%*%t(R[,im])
    CM[Range,]=Qij
    Range<-Range+n.comp
    
    if (im>1)
       { 
        for (jm in (1:(im-1)))
        {
        Xjm<- X[,jm]
        Qij<-t((Xim *Xjm) %*% t(scale2) * X) %*% X - R[,im]%*%t(R[,jm])- R[,jm]%*%t(R[,im])
        CM[Range,]=sqrt(2)*Qij
        Range<-Range+n.comp
        }
       }
    }
    
    V<-rjd.fortran(CM)$V 
     
    B1<-t(V)%*%W
    A<-iW%*%V
    
    keys<-order(colSums(A^2),decreasing=TRUE)
    B<-B1[keys,]
    if (n.comp>1)
    {
    signs.B1<-ifelse(B[,1]<0,-1,1)
    B<-diag(signs.B1)%*% B
    }
    else
    {
    B<-sign(B[1])*B
    }  
    
  if(n.comp == X.cols) A <- solve(B) else A <- t(B) %*% solve(B %*% t(B))
  S <- data.X %*% t(B)
  colnames(S) <- paste("IC.", 1:n.comp, sep="")
  res<-list(A=A, W=B, S=S, Xmu=Col.center)
  class(res) <- "bss"
  return(res)
  }
  
  
  
  
