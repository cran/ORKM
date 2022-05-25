#' Caculate the pardon matrix and the estimator on the RKMeans
RKMeans=function(X,K,V,yita,r,max.iter,truere,method=0){
#' param X is the data matrix
#' param K is the number of cluster  
#' param yita is the regularized parameter
#' param r is the banlance parameter
#' param V is the view of X
#' param max.iter is the max iter
#' param truere is the ture label in data set
#' param method is the caluate the NMI
#' 
#' @return NMI,weight,center,result
#' @export
#'
#' @examples  
#'  yita=0.5;V=2;K=3;r=0.5;max.iter=10;n1=n2=n3=70
#'  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
#'  Xv<-c(X1,X2,X3)
#'  data<-matrix(Xv,n1+n2+n3,2)
#'  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
#'  X<-matrix(data[,1],n1+n2+n3,1) 
#'  truere=data[,2]
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
#'  RKMeans(X=X1,K=K,V=V,yita=yita,r=r,max.iter=max.iter,truere=truere,method=0)

 if (V<=1){
## RKMeans for single-view
rows<- nrow(X) 
cols<- ncol(X)
ukiMatrix <- matrix(0,nrow=rows,ncol=2)
changed=2
iter=0
 M1 <- matrix(0,nrow=K,ncol=cols)
 M11 <- matrix(0,nrow=K,ncol=cols)
 SJS <- as.vector(sample(1:rows,size=K))
  for (k in 1:K) { 
    ukiMatrix[SJS[k],1] <- k
    M1[k,] <- X[SJS[k],] 
    M1 <- matrix(M1,K,cols)
      }
  while(changed>0.00000001){ 
       M11=M1
       if(iter >= max.iter)
        break
       C1<-matrix(0,rows,K)
       U<-matrix(0,rows,K)
       C2<-matrix(0,rows,1)
      
       for (i in 1:rows){
              for(k in 1:K){
                   J<-matrix(X[i,]- M1[k,],nrow=1,ncol=cols)
                   C1[i,k]<-exp(-(J%*%t(J))/yita)
                     }      
        }
       C2<-apply(C1,1,sum)
       U<-C1/C2
            M2<-matrix(0,K,cols)
           M3<-c(0,0,0)
            for (k in 1:K){
             M2[k,]<-t(U[,k])%*%X
             M3<-apply(U,2,sum)
             M1[k,]<-M2[k,]/M3[k]
           }  
for(i in 1:rows){
           MaxD <- 10000
           Ji <- ukiMatrix[i,1]        
             for(k in 1:K){ 
               J1<-matrix(X[i,]- M1[k,],nrow=1,ncol=cols)
               d <- sum((J1%*%t(J1))/U[i,k])
                 if(d < MaxD){
                    MaxD <- d
                    ukiMatrix[i,1] <- k
                    ukiMatrix[i,2] <- d
                  }
             }
         }
        iter=(iter+1) 
 changed=norm((M11-M1),type="1")
  }}
  if (V>1){
X1<-as.matrix(X) 
N1<-nrow(X1)
J1<-ncol(X1) #2708*2708
iter=1
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
    M1[k,] <- X1[SJS[k],] 
    M1 <- matrix(M1,K,J1)
      }    # Random center matrix 
change=1
u<-matrix(0,1,K)
g<-matrix(0,1,K)
value<-matrix(0,1,K)
P1<-matrix(0,1,N1)

while(change>0.1){
M11=M1
if(iter>=max.iter)
break
Dwave<-Alpha*D   
Dwave[which(Dwave==Inf)]=1
Mt1<- t(X1)%*%Dwave%*%U1%*%(ginv(t(U1)%*%Dwave%*%U1))
M1<-t(Mt1)
for ( i in 1:N1){
    for (k in 1:K){
  g[k]=1
  value[k]=Dwave[k,k]*norm((X1[i,]-g%*%M1),type="2")+yita*norm((g%*%t(g)),type="1")
g<-matrix(0,1,K)
}
  k1<-which.min(value)
 P1[i]<-k1
 u[k1]=1
U1[i,]<-u
u<-matrix(0,1,K)
value<-matrix(0,1,K)
}
value5<-norm((X1- U1%*%M1),type="1")
for (j in 1:N1){
D[j,j]<-(1/(2*value5))
}
value6<-diag(value5)
value7<-sum(value6) 
Alpha1<-(r*value7)^(1/(1-r))
change=norm((M11-M1),type="1")
iter=(iter+1)
}
 }
    if(V<=1){
ccc<-c(ukiMatrix[,1])
weight=1
}else{ccc<-c(P1)
weight=Alpha1
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
list(NMI=NMI,weight=weight,center=M1,result=ccc)
)}



