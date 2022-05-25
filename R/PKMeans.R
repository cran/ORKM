#' Caculate the pardon matrix and the estimator on the PKMeans
PKMeans=function(X,K,yitapower,sm,max.m,truere,method=0){
#' param X is the data matrix
#' param K is the number of cluster  
#' param yitapower is the regularized parameter
#' param sm is the banlance parameter
#' param chushi is the initial value for online
#' param max.m is the max iter
#' param truere is the ture label in data set
#' param method is the caluate the NMI
#' 
#' @return center,NMI,result
#' @export
#'
#' @examples   
#'  yitapower=0.5;K=3;sm=0.5;max.m=100;n1=n2=n3=70
#'  X1<-rnorm(n1,20,2);X2<-rnorm(n2,25,1.5);X3<-rnorm(n3,30,2) 
#'  Xv<-c(X1,X2,X3)
#'  data<-matrix(Xv,n1+n2+n3,2)
#'  data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
#'  truere=data[,2]
#'  X11<-matrix(data[,1],n1+n2+n3,1) 
#'  PKMeans(X=X11,K=K,yitapower=yitapower,sm=sm,max.m=max.m,truere=truere,method=0)
 rows <- nrow(X) 
 cols <- ncol(X)
 m = 1
 changed=2
 powerMatrix <- matrix(0,nrow=rows,ncol=2) 
 M <- matrix(0,nrow=K,ncol=cols)
 M1 <- matrix(0,nrow=K,ncol=cols)
 SJS <- as.vector(sample(1:rows,size=K))
   for (k in 1:K) { 
    powerMatrix[SJS[k],1] <- k
    M[k,] <- X[SJS[k],] 
    M <- matrix(M,K,cols)
  }
   while(changed>0.000000001){
       M1=M
       sm=-0.5
         if (m>=max.m)
         break
        dM1<-matrix(0,rows,K)  
        wM<-matrix(0,rows,K)
       for (i in 1:rows){
           for(j in 1:K){
         X1<-matrix(X[i,],nrow=1,ncol=cols)
         dM1[i,j]<-(X1-M[j,])%*%t(X1-M[j,]) 
                 }
                }
       dM11<-matrix(apply(dM1,1,sum),rows,1)
       wM2<-dM1^(2*(sm-1))
       wM2[which(wM2==Inf)]=1
       wM3<-(dM11)^(2*sm)
       wM1<-(wM3)^(1/sm-1)
       xwM1<-as.vector(t(wM1))
            for (g in 1:rows){
      wM[g,]<- xwM1[g]%*%wM2[g,]
            }
             M2<-matrix(0,K,cols)
             M3<-c(0,0,0)
            for (k in 1:K){
     M2[k,]<-(t(wM[,k])%*%X) 
     M3<-apply(wM,2,sum)
    M[k,]<- M2[k,]/M3[k]
            }
 for(i in 1:rows){
           MaxD <- Inf      
             for(k in 1:K){ 
                  J<-matrix(X[i,]- M[k,],nrow=1,ncol=cols)
               d <- sum((J%*%t(J))/wM[i,k])
                 if(d < MaxD){
                    MaxD <- d
                    powerMatrix[i,1] <- k
                    powerMatrix[i,2] <- d
                  }
             }
}
    m=(m+1)
    sm=sm*yitapower
changed=norm((M1-M),type="1")
   }
if(method==0){
ccc<-powerMatrix[,1]
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
list(center=wM,NMI=NMI,result=ccc))}


