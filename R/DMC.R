#' Deep matrix clustering algorithm for multi-view data
DMC=function(X,K,V,r,lamda,truere,max.iter,method=0){
#' param X is the multi-view data matrix
#' param K is the number of clusters in the input data matrix  
#' param lamda is the parameter of the depth matrix in the DMC algorithm
#' param r is the banlance parameter of the algorithm
#' param V is the total views of X
#' param max.iter is the maximum number of iterations of the algorithm
#' param truere is the true label vector for the calculated dataset
#' param method refers to the calculation of the clustering evaluation indicator NMI
#' 
#' @return NMI,Alpha1,center,result
#' @export
#'
#' @examples  
#'  V=2;lamda=0.5;K=3;r=0.5;max.iter=10;n1=n2=n3=70
#'  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
#'  Xv<-c(X1,X2,X3)
#'  data<-matrix(Xv,n1+n2+n3,2)
#'  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
#'  truere=data[,2]
#'  X<-matrix(data[,1],n1+n2+n3,1) 
#'  lamda1<-0.2;lamda2<-0.8
#'  lamda0<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
#'  sol.svd <- svd(lamda0)
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
#'  DMC(X=X1,K=K,V=V,r=r,lamda=lamda,truere=truere,max.iter=max.iter,method=0)

X11<-as.matrix(X) 
X1<-t(X11)
N1<-nrow(X1)
J1<-ncol(X1) 
alpha<-1/V 
M1 <- matrix(0,nrow=N1,ncol=K)
 M11 <- matrix(0,nrow=N1,ncol=K)

 SJS <- as.vector(sample(1:J1,size=K))
for (k in 1:K) { 
    M1[,k] <- X1[,SJS[k]] 
    M1 <- matrix(M1,N1,K)
      }    
change=1
HM<-matrix(0,nrow=K,ncol=J1)
for(i in 1:J1 )  {
 mr=sample(1:K,1,replace=FALSE)
HM[mr,i]=1
} 
HMI<-matrix(0,nrow=N1,ncol=J1) 

HMI<-M1%*%HM   
Mfront<-ginv(t(M1)%*%M1)  
Mmid<-(alpha^r)*t(M1)%*%X1%*%t(HMI) 
Mlast<-ginv((alpha^r)*HMI%*%t(HMI))  
M1<-t(Mfront%*%Mmid%*%Mlast)
Q1<-matrix(0,nrow=K,ncol=J1)
Q1<-(alpha)^r*t(M1)%*%X1    
P1<-matrix(0,nrow=K,ncol=K)
P1<-alpha^r*t(M1)%*%M1  

FM<-matrix(0,nrow=K,ncol=J1)
FM<-HM%*%t(HM)%*%HM   
Q1true<- matrix(0,nrow=K,ncol=J1)
Q1false<- matrix(0,nrow=K,ncol=J1)
for (i in 1:K){
  for (j in 1:J1){
Q1true[i,j]=(abs(Q1[i,j])+Q1[i,j])/2 
Q1false[i,j]= (abs(Q1[i,j])-Q1[i,j])/2 
}
}
P1true<- matrix(0,nrow=K,ncol=K)
P1false<- matrix(0,nrow=K,ncol=K)
for (i in 1:K){
  for (j in 1:K){
P1true[i,j]=(abs(P1[i,j])+P1[i,j])/2 
P1false[i,j]= (abs(P1[i,j])-P1[i,j])/2 
}
}
FMtrue<-matrix(0,nrow=K,ncol=J1)
FMfalse<-matrix(0,nrow=K,ncol=J1)
for (i in 1:K){
  for (j in 1:J1){
FMtrue[i,j]=(abs(FM[i,j])+FM[i,j])/2 
FMfalse[i,j]= (abs(FM[i,j])-FM[i,j])/2 
}
}
HMlast<-matrix(0,nrow=K,ncol=J1) 
HMlast<-sqrt(Q1true+ P1false%*%HM+lamda* FMfalse)/(Q1false+P1true%*%HM+lamda* FMtrue) 
HMlast[is.na(HMlast)]<-0
HM<-HM*HMlast  
HM[is.na(HM)]<-0
thtaM<-norm((X1- M1%*%HM),type="1") 
Alpha1<-(r* thtaM)^(1/(1-r))  



for (j in 1:J1){
 for (k in 1:K){
 if(k==which.max(HM[,j])){
   HM[k,j]=1}else{
HM[k,j]=0}
}
}
g1<-c(0,J1)
for (i in 1: J1){
g1[i]<-which.max(HM[,i])
}
if(method==0){
ccc<-c(g1)
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
list(NMI=NMI,Alpha1=Alpha1,center=M1,result=ccc)
)}

