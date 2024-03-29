\name{RKMeans}
\alias{RKMeans}
\title{
Regularized K-means clustering algorithm for multi-view data
}
\description{
This function improves the regularized K-means clustering (RKMC) algorithm for the multi-view data clustering problem. Specifically, the regularisation term is added to the K-means algorithm to avoid overfitting of the data. Numerical analysis shows that the RKMC algorithm significantly improves the clustering performance compared to other methods. In addition, in order to reveal the structure of real data as realistically as possible, improve the clustering accuracy of high-dimensional data, and balance the weights of each view, the RKMC algorithm assigns a series of learnable weight values to each view, thus reflecting the relationship and compatibility of each view more flexibly.
}
\usage{
RKMeans(X, K, V, yita, r, max.iter, truere, method = 0)
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
  \item{yita}{
is the regularized parameter
}
  \item{r}{
is the banlance parameter
}
  \item{max.iter}{
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
NMI,weight,center,result
}
\author{
Miao Yu
}
\examples{
  library(MASS) 
  library(Matrix)  
  yita=0.5;V=2;K=3;r=0.5;max.iter=10;n1=n2=n3=70
  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
  Xv<-c(X1,X2,X3)
  data<-matrix(Xv,n1+n2+n3,2)
  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
  X<-matrix(data[,1],n1+n2+n3,1) 
  truere=data[,2]
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
  RKMeans(X=X1,K=K,V=V,yita=yita,r=r,max.iter=max.iter,truere=truere,method=0)
}