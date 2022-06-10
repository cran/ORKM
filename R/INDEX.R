#' Caculate the indication on the functions
INDEX = function(vec1, vec2,method=0,mybeta=0) {
#' param vec1: the cluster result by experiment
#' param vec2: the real cluster, stand cluster  
#' param method=0 purity
#' param method=1 precision
#' param method=2 recall
#' param method=3 F-score with para mybeta
#' param method=4 RI
#' 
#' @return accuracy
#' @export
#'
#' @examples 
#' P1<-c(1,1,1,2,3,2,1);truelabel<-c(1,1,1,2,2,2,3)
#' INDEX(P1,truelabel,method=0);INDEX(P1,truelabel,method=2)
    n = length(vec1)
    
    nodes = 1:n
    
    K1= length(unique(vec1))
    K2= length(unique(vec2))
    com1 = c()
    com2 = c()
    for (i in 1:K1) {
       com1[i] = list(nodes[vec1 == i])
    }
for (i in 1:K2) {
       com2[i] = list(nodes[vec2 == i])
    }
    # browser()
    if(method==0){
## purity
    cor_num = 0
    may_id = 1:K2
    for (i in 1:K1) {
       
       temp = sapply(may_id, function(x) {
          return(length(intersect(com1[[i]], com2[[x]])))
        
       })
       cor_num=cor_num+max(temp)   
       # id_del = which.max(temp)
       # cor_num = cor_num + max(temp)
       # may_id = may_id[-id_del]
       
    }
    
    accuracy = cor_num/n
    }
if(method!=0){
## precision
expe_cluster_mat=matrix(0,n,n);
real_cluster_mat=matrix(0,n,n);
for(i in 1:K1){
nodepair_id=combn(com1[[i]],2);
# print(dim(nodepair_id)[2])
node_id=apply(nodepair_id,2,function(x){
                     x=sort(x);
 x[1]+(x[2]-1)*n
                          })
expe_cluster_mat[node_id]=1;
# length(which(expe_cluster_mat==1))
}
expe_cluster_mat[lower.tri(expe_cluster_mat)] =2;
# expe_cluster_mat[lower.tri(expe_cluster_mat)] = t(expe_cluster_mat)[lower.tri(t(expe_cluster_mat))]
    for(i in 1:K2){
nodepair_id=combn(com2[[i]],2);
node_id=apply(nodepair_id,2,function(x){
                     x=sort(x);
 x[1]+(x[2]-1)*n
                          })
real_cluster_mat[node_id]=1;
}
# real_cluster_mat[lower.tri(real_cluster_mat)] = t(real_cluster_mat)[lower.tri(t(real_cluster_mat))]
real_cluster_mat[lower.tri(real_cluster_mat)] =2;
diag(expe_cluster_mat)=2;
diag(real_cluster_mat)=2;
# browser()
exp_p=which(expe_cluster_mat==1)
rel_p=which(real_cluster_mat==1)
exp_n=which(expe_cluster_mat==0)
rel_n=which(real_cluster_mat==0)
# choose(8,2)+choose(5,2)+choose(4,2)
TP=length(intersect(rel_p,exp_p))
# FP=length(exp_p)-TP;
FP=length(intersect(rel_n,exp_p))

## -n, delete the effect by the diag element, which both be 0
TN=length(intersect(rel_n,exp_n))
FN=length(intersect(rel_p,exp_n))
# FN=length(rel_n)-TN;
P=TP/(TP+FP);
R=TP/(TP+FN);
} ## end if method !=0
 

if(method==1){
accuracy=P
}else if (method==2){
accuracy=R
}else if(method==3){
accuracy=(mybeta+1)*P*R/(mybeta*P+R)
}else if(method==4){
accuracy=(TP+TN)/(TP+TN+FP+FN)
}
    return(accuracy)
} 
