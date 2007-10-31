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
    if (n.comp > X.cols) stop("'n.comp' must be smaller than number of colums of X")
    
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
    
    V<-diag(1,n.comp)
    encore<-1
    sweep<-0
    updates<-0
    g<-matrix(0,ncol=2,nrow=nbcm)
    gg<-matrix(0,ncol=2,nrow=2)
    G<-matrix(0,ncol=2,nrow=2)
    ton<-0
    toff<-0
    theta<-0

    if (n.comp>1)
    {
    while (encore==1)
        {
        encore<-0
        for (p in 1:(n.comp-1))
            {
            for (q in (p+1):n.comp)
                {
        
                Ip<-seq(p,n.comp*nbcm,n.comp)
                Iq<-seq(q,n.comp*nbcm,n.comp)
        
                g[,1]<-CM[Ip,p]-CM[Iq,q]
                g[,2]<-CM[Iq,p]+CM[Ip,q]
            
                gg<-t(g)%*%g
                ton<-gg[1,1]-gg[2,2]
                toff<-gg[1,2]+gg[2,1]
                theta<-0.5*atan2(toff, ton+sqrt(ton*ton+toff*toff))
        
                if (abs(theta)>eps)
                    {
                    encore<-1
                    updates<-updates+1
                    cos.theta<-cos(theta)
                    sin.theta<-sin(theta)
                    G<-rbind(c(cos.theta,-sin.theta),c(sin.theta,cos.theta))
                    pair<-c(p,q)
                    V[,pair]<-V[,pair] %*% G
                    CM[,pair]<-CM[,pair] %*% G
                    CM[c(Ip,Iq),]<-rbind(cos.theta*CM[Ip,]+sin.theta*CM[Iq,], -sin.theta*CM[Ip,]+cos.theta*CM[Iq,])
                    }
                }
            }
        if (updates >= maxiter) stop("maxiter reached without convergence")   
        }
    }
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
  res<-list(A=solve(B),W=B,S=data.X %*% t(B),Xmu=Col.center)
  return(res)
  }

