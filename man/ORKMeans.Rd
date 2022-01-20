\name{ORKMeans}
\alias{ORKMeans}
\title{Caculate the estimator on the ORKMeans method}
\description{
Caculate the estimator on the ORKMeans method
}
\usage{
ORKMeans(X, K, V, chushi, yita, gamma, truere, max.iter, method = 0)
}
\arguments{
  \item{X}{
is the data matrix
}
  \item{K}{
is the number of cluster 
}
  \item{V}{
is the view of X
}
  \item{chushi}{
is the initial value for online
}
  \item{yita}{
is the regularized parameter
}
  \item{gamma}{
is the banlance parameter
}
  \item{truere}{
is the ture label in data set
}
  \item{max.iter}{
is the max iter
}
  \item{method}{
is the caluate the NMI
}
}
\value{
mvNM,mvAlpha1,mvonM,mvresult
}

\author{
Miao Yu
}
\examples{
library(MASS) 
library(Matrix)  
  yita=0.5;V=2;chushi=100;K=3;gamma=0.5;max.iter=10;n1=n2=n3=70
  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
  Xv<-c(X1,X2,X3)
  data<-matrix(Xv,n1+n2+n3,2)
  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
  truere=data[,2]
  X<-matrix(data[,1],n1+n2+n3,1) 
  lamda1<-0.2;lamda2<-0.8
  lamda<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
  sol.svd <- svd(lamda)
  U1<-sol.svd$u
  D1<-sol.svd$d
  V1<-sol.svd$v
  C1<-t(U1)%*%t(X)
  Y1<-C1/D1
  view<-V1%*%Y1
  view1<-matrix(view[1,])
  view2<-matrix(view[2,])
  X1<-matrix(view1,n1+n2+n3,1)
  X2<-matrix(view2,n1+n2+n3,1)
  ORKMeans(X=X1,K=K,V=V,chushi=chushi,yita=yita,gamma=gamma,
max.iter=max.iter,truere=truere,method=0)
}
