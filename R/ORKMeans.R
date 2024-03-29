#' Online regularized K-means clustering algorithm for online multi-view data
ORKMeans=function(X,K,V,chushi,r,yita,gamma,alpha,epsilon,truere,max.iter,method=0){
#' param X is the online mulit-view data matrix
#' param K is the number of clusters in the input data matrix 
#' param yita is the regularization parameter of the algorithm
#' param r is the banlance parameter of the algorithm
#' param gamma is the step size of the algorithm
#' param alpha is the weight of the calculated view
#' param V is the total views of X
#' param chushi is the initial value for online algorithm
#' param epsilon is the algorithm stopping threshold
#' param max.iter is the maximum number of iterations of the algorithm
#' param truere is the true label vector for the calculated dataset
#' param method refers to the calculation of the clustering evaluation indicator NMI
#' 
#' @return NMI,weight,center,result
#' @export
#'
#' @examples   
#'  yita=0.5;V=2;chushi=100;K=3;r=0.5;max.iter=10;n1=n2=n3=70;gamma=0.1;alpha=0.98;epsilon=1
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
#'  ORKMeans(X=X1,K=K,V=V,chushi=chushi,r=r,yita=yita,gamma=gamma,alpha=alpha,epsilon=epsilon,max.iter=max.iter,truere=truere,method=0)
 if (V<=1){
## KMeans for single-view
N<- nrow(X) 
J<- ncol(X)
oX<-matrix(X[1:chushi,],chushi,J)
oN<- nrow(oX) 
oJ<- ncol(oX)
oukiMatrix <- matrix(0,nrow=oN,ncol=2)
changed=2
iter=0


 oM <- matrix(0,nrow=K,ncol=oJ)
 oM1 <- matrix(0,nrow=K,ncol=oJ)
 SJS <- as.vector(sample(1:oN,size=K))
  for (k in 1:K) { 
    oukiMatrix[SJS[k],1] <- k
    oM[k,] <- X[SJS[k],] 
    oM <- matrix(oM,K,oJ)
      }
  while(changed>0.00000001){ 
       oM1=oM
       if(iter >= max.iter)
        break
       oC1<-matrix(0,oN,K)
       oU<-matrix(0,oN,K)
       oC2<-matrix(0,oN,1)    
       for (i in 1:oN){
              for(k in 1:K){
                   ooJ<-matrix(oX[i,]- oM[k,],nrow=1,ncol=oJ)
                   oC1[i,k]<-exp(-(ooJ%*%t(ooJ))/yita)
                     }      
        }
       oC2<-apply(oC1,1,sum)
       oU<-oC1/oC2
            oM2<-matrix(0,K,oJ)
           oM3<-c(0,0,0)
            for (k in 1:K){
             oM2[k,]<-t(oU[,k])%*%oX
             oM3<-apply(oU,2,sum)
             oM[k,]<-oM2[k,]/oM3[k]
           }  
for(i in 1:oN){
           MaxD <- 10000
           Ji <- oukiMatrix[i,1]        
             for(k in 1:K){ 
               J1<-matrix(oX[i,]- oM[k,],nrow=1,ncol=oJ)
               d <- sum((J1%*%t(J1))/oU[i,k])
                 if(d < MaxD){
                    MaxD <- d
                    oukiMatrix[i,1] <- k
                    oukiMatrix[i,2] <- d
                  }
             }
         }
        iter=(iter+1) 
 changed=norm((oM1-oM),type="1")
  }

onM<-oM
  onC1<-matrix(0,N,K)
  onU<-matrix(0,N,K)
  onC2<-c(0,N,1)
P2<-matrix(0,1,N)
ling<-matrix(0,N-chushi,K)
onU<-rbind(oU,ling)

for (i in (chushi+1):N){
    for(k in 1:K){ 
 oo<-matrix(X[i,]-onM[k,],nrow=1,ncol=J) 
          onC1[i,k]<-onU[i-1,k]%*%exp(-(oo%*%t(oo))/yita)
 onC2[i]<-sum(onC1[i,])
     onU[i,]<-onC1[i,]/onC2[i]
      }
onU[is.na(onU)]<-0
       onM2<-matrix(0,K,J)
  onk1<-which.max(onU[i,])
 P2[i]<-onk1
     onM[onk1,]<-(onM[onk1,]+X[i,])/2
}

oP2<-matrix(0,1,N)
for (i in 1:chushi){
  ok1<-which.max(oU[i,])
 P2[i]<-ok1
}

}

######### KMeans for multi-view
  if (V>1){

N1<-nrow(X)
J1<-ncol(X)
cX<-matrix(X[1:chushi,],chushi,J1)
cN1<-nrow(cX)
cJ1<-ncol(cX)
iter=0  
alpha<-1/V  
D<-diag(cN1)   

cU1<-matrix(0,cN1,K) 
cU1<-matrix(0,nrow=cN1,ncol=K)
for(i in 1:cN1)  {
 mr=sample(1:K,1,replace=FALSE)
cU1[i,mr]=1
}  
M1 <- matrix(0,nrow=K,ncol=cJ1)
M11 <- matrix(0,nrow=K,ncol=cJ1)

 SJS <- as.vector(sample(1:cN1,size=K))

for (k in 1:K){ 
    M1[k,] <- cX[SJS[k],] 
    M1 <- matrix(M1,K,cJ1)
      }     
change=1
cu<-matrix(0,1,K)
cg<-matrix(0,1,K)
cvalue<-matrix(0,1,K)
cP1<-matrix(0,1,cN1)

while(change>0.1){
M11=M1
if(iter>=max.iter)
break
Dwave<-alpha*D  
 
Dwave[which(Dwave==Inf)]=1
Mt1<- t(cX)%*%Dwave%*%cU1%*%(ginv(t(cU1)%*%Dwave%*%cU1))
M1<-t(Mt1)
for ( i in 1:cN1){
    for (k in 1:K){
  cg[k]=1
  cvalue[k]=Dwave[k,k]*norm((cX[i,]-cg%*%M1),type="2")+yita*norm((cg%*%t(cg)),type="1")
cg<-matrix(0,1,K)
}
  ck1<-which.min(cvalue)
 cP1[i]<-ck1
 cu[ck1]=1
cU1[i,]<-cu
cu<-matrix(0,1,K)
cvalue<-matrix(0,1,K)
}
value5<-norm((cX- cU1%*%M1),type="1")
for (j in 1:cN1){
D[j,j]<-(1/(2*value5))
}
value6<-diag(value5)
value7<-sum(value6) 
Alpha1<-(gamma*value7)^(1/(1-gamma))
change=norm((M11-M1),type="1")
iter=(iter+1)
}

alpha<-alpha
oM1<-matrix(0,cJ1,K)
oM1<-t(M1)
A<-matrix(0,N1,K) 
ochange=1
gamma=gamma
r=r
oU1<-matrix(0,nrow=N1,ncol=K)
onU1<-matrix(0,nrow=N1,ncol=K)
pU1<-matrix(0,nrow=N1,ncol=K)
pu<-c(rep(0,K))
P1<-c(rep(0,N1))

D<-diag(N1) 
Dwave<-alpha*D 
for (i in (chushi+1):N1){
ling<-matrix(0,N1-chushi,K)
onU1<-rbind(cU1,ling)
g = 1
A[1:i,]<-X[1:i,]%*%oM1
epsilon <- epsilon

dJ <- function(onU1){
  onU1[which(onU1==-Inf)]=0 
  return((2* Dwave[1:i,1:i]%*% onU1[1:i,]%*%t(oM1) %*%oM1+2*yita*onU1[1:i,]-2*Dwave[1:i,1:i]%*% A[1:i,]))
}
J <- function(onU1){
  onUU<-onU1[1:i,]%*%t(onU1[1:i,])
  onUU[which(onUU==Inf)]=1
  return(alpha*(norm((X[1:i,]-onU1[1:i,]%*%t(oM1)),type="1")+yita*sum(diag(onUU))))
}

while(TRUE){
  gradient = dJ(onU1)
  last_onU1 = onU1
  onU1[which(onU1==-Inf)]=0 
  onU1[1:i,] = onU1[1:i,] - gamma * gradient
  onU1[which(onU1==-Inf)]=0
  g <- g+1
  
  if (abs(J(onU1) - J(last_onU1)) < epsilon){
    break
  }

}
}

for (i in (chushi+1):N1){
k1<-which.max(onU1[i,])
P1[i]<-k1 
pu[k1]<-1
pU1[i,]<-pu
pu<-c(rep(0,K))
onM1= t(X[1:i,])%*%Dwave[1:i,1:i]%*%pU1[1:i,]%*%ginv(t(pU1[1:i,])%*%Dwave[1:i,1:i]%*%pU1[1:i,])

}
for (i in 1:chushi){
pU1[i,]<-cU1[i,]
P1[i]<-cP1[i]
}
value<-c(rep(0,N1))
for (i in chushi:N1){
value [i]<-r*norm((X[1:i,]- pU1[1:i,]%*%t(onM1)),type="1")  
}
value1<-sum(value^(1/(1-r)))  
Alpha1<-(r* value[i]^(1/(1-r)) /value1)
}
 if(V<=1){
ccc<-c(P2)
reM=onM
alpha=1
}else{ccc<-c(P1)
reM=onM1
alpha=Alpha1
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
list(NMI=NMI,weight=alpha,center=reM,result=ccc)
)

}
