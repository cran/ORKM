\name{ORKM-package}
\alias{ORKM-package}
\alias{ORKM}
\docType{package}
\title{
\packageTitle{ORKM}
}
\description{
\packageDescription{ORKM}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ORKM}
\packageIndices{ORKM}
You can use this package for online multi-view clustering, the dataset and real labels are also provided in the package. 
}
\author{
\packageAuthor{ORKM}

Maintainer: \packageMaintainer{ORKM}
}
\references{
Guangbao Guo, Miao Yu, Guoqi Qian, (2023), Orkm: Online Regularized k-Means Clustering for Online Multi-View Data.
}

\keyword{online multi-view clustering}
\seealso{
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4484209
}
\examples{
library(MASS) 
library(Matrix)  
  yita=0.5;V=2;chushi=100;K=3;r=0.5;max.iter=10;n1=n2=n3=70;gamma=0.1;alpha=0.98;epsilon=1
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
  ORKMeans(X=X1,K=K,V=V,r=r,chushi=chushi,yita=yita,gamma=gamma,epsilon=epsilon,
max.iter=max.iter,truere=truere,method=0)
}
