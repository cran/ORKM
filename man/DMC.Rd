\name{DMC}
\alias{DMC}
\title{Deep matrix clustering algorithm for multi-view data}
\description{
This algorithm decomposes the multi-view data matrix into representative subspaces layer by layer, and generates a cluster at each layer. To enhance the diversity between the generated clusters, new redundant quantifiers arising from the proximity between samples in these subspaces are minimised. An iterative optimisation process is further introduced to simultaneously seek multiple clusters with quality and diversity.
}
\usage{
DMC(X, K, V, r, lamda, truere, max.iter, method = 0)
}
\arguments{
  \item{X}{data matrix}

  \item{K}{number of cluster}

  \item{V}{number of view}

  \item{r}{first banlance parameter}

  \item{lamda}{second balance parameter}

  \item{truere}{true cluster result}

  \item{max.iter}{max iter}

  \item{method}{caculate the index of NMI}
}
\value{
NMI,Alpha1,center,result
}
\author{
Miao Yu
}
\examples{
 library(MASS)   
 V=2;lamda=0.5;K=3;r=0.5;max.iter=10;n1=n2=n3=70
 X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
 Xv<-c(X1,X2,X3)
 data<-matrix(Xv,n1+n2+n3,2)
 data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
 truere=data[,2]
 X<-matrix(data[,1],n1+n2+n3,1) 
 lamda1<-0.2;lamda2<-0.8
 lamda0<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
 sol.svd <- svd(lamda0)
 U1<-sol.svd$u
 D1<-sol.svd$d
 V1<-sol.svd$v
 C1<-t(U1)\%*\%t(X)
 Y1<-C1/D1
 view<-V1\%*\%Y1
 view1<-matrix(view[1,])
 view2<-matrix(view[2,])
 X1<-matrix(view1,n1+n2+n3,1)
 X2<-matrix(view2,n1+n2+n3,1)
 DMC(X=X1,K=K,V=V,lamda=lamda,r=r,max.iter=max.iter,truere=truere,method=0)
}
