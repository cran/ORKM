#' Caculate the pardon matrix and the estimator on the online RKMeans
ORKMeans=function(X,K,V,chushi,yita,gamma,truere,max.iter,method=0){
#' param X is the data matrix
#' param K is the number of cluster  
#' param yita is the regularized parameter
#' param gamma is the banlance parameter
#' param V is the view of X
#' param chushi is the initial value for online
#' param max.iter is the max iter
#' param truere is the ture label in data set
#' param method is the caluate the NMI
#' 
#' @return mvNM,mvAlpha1,mvonM,mvresult
#' @export
#'
#' @examples   
#'  yita=0.5;V=2;chushi=100;K=3;gamma=0.5;max.iter=10;n1=n2=n3=70
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
#'  ORKMeans(X=X1,K=K,V=V,chushi=chushi,yita=yita,gamma=gamma,max.iter=max.iter,truere=truere,method=0)
 if (V<=1){
## KMeans for single-view
N<- nrow(X) 
J<- ncol(X)
oX<-matrix(X[1:chushi,],chushi,1)
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
for (i in (chushi+1):N){
    for(k in 1:K){ 
 oo<-matrix(X[i,]-onM[k,],nrow=1,ncol=J) 
          onC1[i,k]<-exp(-(oo%*%t(oo))/yita)
      }
 onC2[i]<-sum(onC1[i,])
     onU[i,]<-onC1[i,]/onC2[i]
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
      }    # Random center matrix 
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

alpha<-0.98
cM1<-matrix(0,chushi,K)
cM1<-t(M1)
ochange=1
tU1<-matrix(0,nrow=N1,ncol=K)
for(i in 1:N1 )  {
 mr=sample(1:K,1,replace=FALSE)
tU1[i,mr]=1
}  
u<-matrix(0,1,K)
g<-matrix(0,1,K)
value<-matrix(0,1,K)
P1<-matrix(0,1,N1)
Mt1<-matrix(0,N1,K)
for (i in (chushi+1):N1){
D<-diag(N1)
Dwave<-alpha*D  
Mt1<- t(matrix(X[i,],1,J1))%*%matrix(Dwave[i,i],1,1)%*%tU1[i,]%*%ginv(t(tU1)%*%Dwave%*%tU1)
onM<-t(Mt1) #online update clustering center matrix 
    for (k in 1:K){
  g[k]=1
  value[k]=alpha*norm((t(X[i,])-g%*%onM),type="2")+yita*norm((g%*%t(g)),type="1")
g<-matrix(0,1,K)
}
  k1<-which.min(value)
 P1[i]<-k1
 u[k1]=1
tU1[i,]<-u
u<-matrix(0,1,K)
value<-matrix(0,1,K)
value5<-c(0,1,J1-(chushi+1))
value5[i]<-norm((t(X[i,])- tU1[i,]%*%onM),type="1")
D[i,i]<-(1/(2*value5[i]))
value6<-diag(value5[i])
value7<-sum(value6) 
Alpha1<-(gamma*value7)^(1/(1-gamma))
}
for (i in 1:chushi){
 P1[i]<-cP1[i]
}
}
 if(V<=1){
ccc<-c(P2)
}else{ccc<-c(P1)
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
if(V<=1){
 return(
list(svNMI=NMI,svonM=onM,svresult=ccc)
)
}else {
return(
list(mvNMI=NMI,mvAlpha1=Alpha1,mvonM=onM,mvresult=ccc)
)
}}



