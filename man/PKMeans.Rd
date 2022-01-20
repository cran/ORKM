\name{PKMeans}
\alias{PKMeans}
\title{
Caculate the pardon matrix and the estimator on the PKMeans
}
\description{
Caculate the pardon matrix and the estimator on the PKMeans
}
\usage{
PKMeans(X, K, yitapower, sm, max.m, truere, method = 0)
}
\arguments{
  \item{X}{
is the data matrix
}
  \item{K}{
is the number of cluster 
}
  \item{yitapower}{
is the regularized parameter
}
  \item{sm}{
is the banlance parameter
}
  \item{max.m}{
is the max iter
}
  \item{truere}{
is the ture label in data set
}
  \item{method}{
is the caluate the NMI
}
}
\value{
wM,wNMI
}
\author{
Miao Yu
}
\examples{
library(MASS)   
  yitapower=0.5;K=3;sm=0.5;max.m=100;n1=n2=n3=70
  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
  Xv<-c(X1,X2,X3)
  data<-matrix(Xv,n1+n2+n3,2)
  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
  truere=data[,2]
  X11<-matrix(data[,1],n1+n2+n3,1) 
  PKMeans(X=X11,K=K,yitapower=yitapower,sm=sm,max.m=max.m,truere=truere,method=0)
}
