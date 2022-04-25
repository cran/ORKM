\name{OGD}
\alias{OGD}
\title{Caculate the result and the estimator on the OGD}

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
Caculate the pardon matrix and the estimator on the OGD
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

