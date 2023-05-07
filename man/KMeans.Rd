\name{KMeans}
\alias{KMeans}
\title{K-means clustering algorithm for multi/single view data}
\description{
The K-means clustering algorithm is a common clustering algorithm that divides a data set into K clusters, with each cluster represented using the mean of all samples within the cluster, referring to that mean as the j-cluster centre. The algorithm is unsupervised learning, where the categories are not known in advance and similar objects are automatically grouped into the same cluster. The K-means algorithm achieves clustering by calculating the distance between each point and the centre of mass of different clusters and assigning it to the nearest cluster. The algorithm is simple and easy to implement, but is susceptible to the initial centre of mass, the possibility of empty clusters, and the possibility of convergence to local minima. Clustering applications can be used to discover different groups of users, allowing for tasks such as precision marketing, document segmentation, finding people in the same circle in social networks, and handling anomalous data.
}
\usage{
KMeans(X, K, V, r, max.iter, truere, method = 0)
}
\arguments{  
  \item{X}{data matrix}

  \item{K}{number of cluster}

  \item{V}{number of view}

  \item{r}{balance parameter}

  \item{truere}{true cluster result}

  \item{max.iter}{max iter}

  \item{method}{caculate the index of NMI}

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
  V=2;K=3;r=0.5;max.iter=10;n1=n2=n3=70
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
  KMeans(X=X1,K=K,V=V,r=r,max.iter=max.iter,truere=truere,method=0)
}