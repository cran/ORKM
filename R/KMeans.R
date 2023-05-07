#' K-means clustering algorithm for multi/single view data
KMeans=function(X,K,V,r,max.iter,truere,method=0){
#' param X is the single/multi-view data matrix
#' param K is the number of clusters in the input data matrix  
#' param V is the total views of X
#' param r is the banlance parameter of the algorithm
#' param max.iter is the maximum number of iterations of the algorithm
#' param truere is the true label vector for the calculated dataset
#' param method refers to the calculation of the clustering evaluation indicator NMI
#' 
#' @return NMI,weight,center,result
#' @export
#'
#' @examples  
#'  V=2;K=3;r=0.5;max.iter=10;n1=n2=n3=70
#'  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
#'  Xv<-c(X1,X2,X3)
#'  data<-matrix(Xv,n1+n2+n3,2)
#'  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
#'  truere=data[,2]
#'  X<-matrix(data[,1],n1+n2+n3,1) 
#'  lamda1<-0.2;lamda2<-0.8
#'  lamda<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
#'  sol.svd <- svd(lamda)
#'  U1<-sol.svd$u
#'  D1<-sol.svd$d
#'  V1<-sol.svd$v
#'  C1<-t(U1)%*%t(X)
#'  Y1<-C1/D1
#'  view<-V1%*%Y1
#'  view1<-matrix(view[1,])
#'  view2<-matrix(view[2,])
#'  X1<-matrix(view1,n1+n2+n3,1)
#'  X2<-matrix(view2,n1+n2+n3,1)
#'  KMeans(X=X1,K=K,V=V,r=r,max.iter=max.iter,truere=truere,method=0)
  if (V<=1){
## KMeans for single-view
     rows <- nrow(X) 
     cols <- ncol(X)
     within <- matrix(0,nrow=K,ncol=1) 
     iter = 0
     max.iter=max.iter
     indexMatrix <- matrix(0,nrow=rows,ncol=2) 
     M <- matrix(0,nrow=K,ncol=cols)
     SJS <- as.vector(sample(1:rows,size=K))
  for (i in 1:K) { 
    indexMatrix[SJS[i],1] <- i
    M[i,] <- X[SJS[i],] 
    M <- matrix(M,K,cols)
      }
   changed=TRUE
   while(changed){ 
      if(iter >= max.iter)
        break
      changed=FALSE
      for(i in 1:rows){
           MaxD <- 10000
           Ji <- indexMatrix[i,1]
           for(j in 1:K){ 
             d <- (sum((X[i,]- M[j,])^2))^0.5
             if(d < MaxD){
               MaxD <- d
               indexMatrix[i,1] <- j
               indexMatrix[i,2] <- d
          }
       }
     if(Ji!=indexMatrix[i,1])
        changed=TRUE
  }
   for( m in 1:K){
     clusterMatrix <- X[indexMatrix[,1]==m,] 
     clusterMatrix <- as.matrix(clusterMatrix)
     if(nrow(clusterMatrix)>0){
       M[m,] <- colMeans(clusterMatrix) 
      }
       else{
         M[m,] <- M[m,]
       }
    }
    iter=(iter+1) 
} }
## KMeans for multi-views
  if (V>1){
  N1=nrow(X)
  J1=ncol(X)
  iter=0
  changed=2
  Alpha<-1/V  
  D<-diag(N1)    

U1<-matrix(0,nrow=N1,ncol=K)
for(i in 1:N1 )  {
 mr=sample(1:K,1,replace=FALSE)
U1[i,mr]=1
}  
M1 <- matrix(0,nrow=K,ncol=J1)
M11 <- matrix(0,nrow=K,ncol=J1)

 SJS <- as.vector(sample(1:N1,size=K))

for (k in 1:K) { 
    M1[k,] <- X[SJS[k],] 
    M1 <- matrix(M1,K,J1)
      }     
change=1
u<-matrix(0,1,K)
g<-matrix(0,1,K)
value<-matrix(0,1,K)
P1<-matrix(0,1,N1)
system.time(
while(change>0.1){
M11=M1
if(iter>=max.iter)
break
Dwave<-Alpha*D   
Dwave[which(Dwave==Inf)]=1
Mt1<- t(X)%*%Dwave%*%U1%*%(ginv(t(U1)%*%Dwave%*%U1))
M1<-t(Mt1)
for ( i in 1:N1){
    for (k in 1:K){
  g[k]=1
  value[k]=Dwave[k,k]*norm((X[i,]-g%*%M1),type="2")
g<-matrix(0,1,K)
}
  k1<-which.min(value)
 P1[i]<-k1
 u[k1]=1
U1[i,]<-u
u<-matrix(0,1,K)
value<-matrix(0,1,K)
}
value5<-norm((X- U1%*%M1),type="1")
for (j in 1:N1){
D[j,j]<-(1/(2*value5))
}
value6<-diag(value5)
value7<-sum(value6)
Alpha1<-(r*value7)^(1/(1-r))
change=norm((M11-M1),type="1")
iter=(iter+1)
}
) }
    if(V<=1){
ccc<-c(indexMatrix[,1])
center=M
Alpha=1
}else{ccc<-c(P1)
center=M1
Alpha=Alpha1
}
if(method==0){
kmfrequency<-as.data.frame(table(ccc))  
kf1<-kmfrequency$Freq/length(ccc) 
H_indexre<-(-sum(kf1*log(kf1)))
tfrequency<-as.data.frame(table(truere)) 
kf2<-tfrequency$Freq/length(truere) 
H_truere<-(-sum(kf2*log(kf2)))
cfrequency<-as.data.frame(table(paste(ccc,truere))) 
kf3<-cfrequency$Freq/length(paste(ccc,truere)) 
H_paste<-(-sum(kf3*log(kf3)))
MI<-H_indexre+H_truere- H_paste
NMI<-MI/sqrt(H_indexre* H_truere)
} 
 return(
list(NMI=NMI,weight=Alpha,center=center,result=ccc)
)}
