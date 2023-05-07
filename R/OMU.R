#' Caculate the pardon matrix and the estimator on the OMU
OMU=function(X,K,V,chushi,yita,r,max.iter,epsilon,truere,method=0){
#' param X is the online multi-view data matrix
#' param K is the number of clusters in the input data matrix
#' param V is the total views of X  
#' param chushi is the initial value for online algorithm
#' param yita is the regularization parameter of the algorithm
#' param r is the banlance parameter of the algorithm
#' param max.iter is the maximum number of iterations of the algorithm
#' param epsilon is the algorithm stopping threshold
#' param truere is the true label vector for the calculated dataset
#' param method refers to the calculation of the clustering evaluation indicator NMI
#' 
#' @return NMI,result,M
#' @export
#'
#' @examples   
#'  yita=0.5;V=2;chushi=100;K=3;r=0.5;max.iter=10;n1=n2=n3=70;epsilon=1
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
#'  OMU(X=X1,K=K,V=V,chushi=chushi,yita=yita,r=r,max.iter=max.iter,epsilon=epsilon,truere=truere,method=0)
X1<-as.matrix(X)
N1<-nrow(X1)
J1<-ncol(X1)
cX1<-matrix(X1[1:chushi,],chushi,J1) 

cN1<-nrow(cX1)
cJ1<-ncol(cX1)
iter=0      
alpha<-1/V  
max.iter<-max.iter 
D<-diag(cN1)    
cU1<-matrix(0,cN1,K) 
cU1<-matrix(0,nrow=cN1,ncol=K)
for(i in 1:cN1)  {
 mr=sample(1:K,1,replace=FALSE)
cU1[i,mr]=1
}  
M1 <- matrix(0,nrow=cJ1,ncol=K)
M11 <- matrix(0,nrow=cJ1,ncol=K)
 SJS <- as.vector(sample(1:cN1,size=K))
for (k in 1:K){ 
    M1[,k] <- cX1[SJS[k],] 
    M1 <- matrix(M1,cJ1,K)
      }   
change=1
cu<-matrix(0,1,K)
cg<-matrix(0,1,K)
cvalue<-matrix(0,1,K)
cP1<-matrix(0,1,cN1)
system.time(
while(change>0.1){
M11=M1
if(iter>=max.iter)
break
Dwave<-alpha*D  
 
Dwave[which(Dwave==Inf)]=1
M1<- t(cX1)%*%Dwave%*%cU1%*%(ginv(t(cU1)%*%Dwave%*%cU1))

for ( i in 1:cN1){
    for (k in 1:K){
  cg[k]=1
  cvalue[k]=Dwave[k,k]*norm((cX1[i,]-cg%*%t(M1)),type="2")+yita*norm((cg%*%t(cg)),type="1")
cg<-matrix(0,1,K)
}
  ck1<-which.min(cvalue)
 cP1[i]<-ck1
 cu[ck1]=1
cU1[i,]<-cu
cu<-matrix(0,1,K)
cvalue<-matrix(0,1,K)
}
value5<-norm((cX1- cU1%*%t(M1)),type="1")
for (j in 1:cN1){
D[j,j]<-(1/(2*value5))
}
value6<-diag(value5)
value7<-sum(value6) 
Alpha1<-(r*value7)^(1/(1-r))
change=norm((M11-M1),type="1")
iter=(iter+1)
}
)


oM1<-matrix(0,cJ1,K)
oM1<-M1
onU1<-matrix(0,nrow=N1,ncol=K)
onU1<-matrix(runif(K,0,1),N1,K)

pU1<-matrix(0,nrow=N1,ncol=K)
pu<-c(rep(0,K))
P2<-c(rep(0,N1))
g=0
cheng1<-matrix(0,N1,K)
cheng2<-matrix(0,N1,K)
cheng3<-matrix(0,J1,K)
cheng4<-matrix(0,J1,K)


cheng1=(X1%*%oM1)
cheng2=(onU1%*%t(oM1)%*%oM1)
onU1=onU1*(cheng1/cheng2)

cheng3=t(X1)%*%onU1
cheng4=oM1%*%t(onU1)%*%onU1
oM1=oM1*(cheng3/cheng4)



for (i in (chushi+1):N1){
k1<-which.max(onU1[i,])
P2[i]<-k1 
pu[k1]<-1
pU1[i,]<-pu
pu<-c(rep(0,K))
}
for (i in 1:chushi){
pU1[i,]<-cU1[i,]
P2[i]<-cP1[i]
}
if(method==0){
ccc<-c(P2)
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
list(NMI=NMI,M=oM1,result=P2)
)}


