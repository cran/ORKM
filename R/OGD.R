#' Caculate the pardon matrix and the estimator on the OGD
OGD=function(X,K,gamma,max.m,chushi,yita,epsilon,truere,method=0){
#' param X is the data matrix
#' param K is the number of cluster  
#' param yita is the regularized parameter
#' param gamma is the step size
#' param chushi is the initial value
#' param epsilon is the epsilon of OGD
#' param max.m is the max iter
#' param truere is the ture label in data set
#' param method is the caluate the NMI
#' 
#' @return result,NMI,M
#' @export
#'
#' @examples  
#'  yita=0.5;V=2;K=3;chushi=100;epsilon=1;gamma=0.1;max.m=10;n1=n2=n3=70
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
#'  OGD(X=X1,K=K,gamma=gamma,max.m=max.m,chushi=chushi,yita=yita,epsilon=epsilon,truere=truere,method=0)

changed=1
N<- nrow(X) 
J<- ncol(X)
chushi=chushi
oX<-matrix(X[1:chushi,],chushi,J)
oN<- nrow(oX) 
oJ<- ncol(oX)
K=K
oukiMatrix <- matrix(0,nrow=oN,ncol=2)
iter=0
yita<-yita
max.iter<-max.m

 oM <- matrix(0,nrow=K,ncol=oJ)
 oM1 <- matrix(0,nrow=K,ncol=oJ)
 set.seed(123)
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

pu<-c(rep(0,K))
P2<-c(rep(0,N))
pU<-matrix(0,nrow=N,ncol=K)
onM<-oM
onU<-matrix(0,nrow=N,ncol=K)
pU1<-matrix(0,nrow=N,ncol=K)
onU=oU
ling<-matrix(0,N-chushi,K)
onU<-rbind(oU,ling)
SJS<-runif(N, min = 0, max = 1)
 for (b in (chushi+1):N){
   onU[b,1]<-SJS[b]
   onU[b,2]<-SJS[N-b+1]
}
DJ<-matrix(0,N,K)
J1m<-matrix(0,N,K)

for (i in (chushi):N){
 for (k in 1:K){

dJ <- function(onU){
  return(norm(matrix(X[i,]-onM[k,],1,J),type="1")+yita*(log(onU[i,k])+1))
  }
 DJ[i,k] <- dJ(onU)
 DJ[which(DJ==-Inf)]=0
J1 <- function(onU){
  return(sum(onU[i,k]*norm(matrix(X[i,]-onM[k,],1,J),type="1")+yita*onU[i,k]*log(onU[i,k])))
}
J1m[i,k]<-J1(onU)
J1m[is.na(J1m)] <- 0
g=0
while(TRUE){
  gradient =  DJ[i,k]
  last_onU = onU
  onU[i,k] = onU[i-1,k] - gamma * gradient
  g <- g+1  
  if (abs(J1m[i,k] - J1m[i-1,k]) < epsilon){
    break
  }
}
}
}

for (i in 1:N){
k1<-which.max(onU[i,])
P2[i]<-k1 
pu[k1]<-1
pU1[i,]<-pu
pu<-c(rep(0,K))
}
if(method==0){
ccc<-P2
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
list(result=P2,NMI=NMI,M=onM)
)}

