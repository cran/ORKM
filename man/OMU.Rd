\name{OMU}
\alias{OMU}
\title{Online multiplicative update algorithm for online multi-view data clustering}
\usage{
OMU(X,K,V,chushi,yita,r,max.iter,epsilon,truere,method=0)
}
\arguments{  
  \item{X}{data matrix}

  \item{K}{number of cluster}

  \item{V}{number of view}

  \item{chushi}{the initial value}

  \item{yita}{the regularized parameter}

  \item{r}{banlance parameter}

  \item{max.iter}{max iter}

  \item{epsilon}{epsilon}

  \item{truere}{true cluster result}

  \item{method}{caculate the index of NMI}

}

\value{
NMI,result,M
}
\description{
This algorithm integrates the multiplicative normalization factor as an additional term in the original additivity update rule, which usually has approximately opposite direction. Thus, the improved iteration rule can be easily converted to a multiplicative version. After each iteration After each iteration, non-negativity is maintained.
}
\examples{
 yita=0.5;V=2;chushi=100;K=3;r=0.5;max.iter=10;n1=n2=n3=70;epsilon=1
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
 C1<-t(U1)\%*\%t(X)
 Y1<-C1/D1
 view<-V1\%*\%Y1
 view1<-matrix(view[1,])
 view2<-matrix(view[2,])
 X1<-matrix(view1,n1+n2+n3,1)
 X2<-matrix(view2,n1+n2+n3,1)
 OMU(X=X1,K=K,V=V,chushi=chushi,yita=yita,r=r,max.iter=max.iter,
epsilon=epsilon,truere=truere,method=0)
}
