\name{OGD}
\alias{OGD}
\title{Online gradient descent algorithm for online single-view data clustering}

\usage{
OGD(X, K, gamma, max.m, chushi, yita, epsilon, truere, method = 0)
}
\arguments{  
  \item{X}{data matrix}

  \item{K}{number of cluster}

  \item{gamma}{step size}

  \item{yita}{the regularized parameter}

  \item{truere}{true cluster result}

  \item{max.m}{max iter}

  \item{epsilon}{epsilon}

  \item{chushi}{the initial value}

  \item{method}{caculate the index of NMI}

}
\value{
result,NMI,M
}
\description{
Online gradient descent is an optimisation algorithm in machine learning for when the amount of data is too large to process all the data at the same time. In this algorithm, the model parameters are updated based on a single training sample, rather than using the entire training set. The direction of each update is determined by the direction of the gradient of the current sample, and the local or global extremes of the gradient descent algorithm depend on the order of the sampled samples. Compared to Batch Gradient Descent (BGD) algorithm, online gradient descent algorithms can process data streams and update the model as they process the data, and are therefore more efficient for large-scale data. However, online gradient descent algorithm should only be used if the data stream is continuously present and updated.
}
\author{
Miao Yu
}
\examples{yita=0.5;V=2;K=3;chushi=100;epsilon=1;gamma=0.1;max.m=10;n1=n2=n3=70
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
 OGD(X=X1,K=K,gamma=gamma,max.m=max.m,chushi=chushi,
yita=yita,epsilon=epsilon,truere=truere,method=0)
}

