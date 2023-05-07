\name{cora_view1}
\alias{cora_view1}
\docType{data}
\title{
The first view of Cora data set. 
}
\description{
This data matrix is the first view of the multi-view data set called Cora, the keyword view. The Cora data set is a multi-view data set of machine learning papers with 4 views, a sample size of nearly 3000 and a number of features of 1500, with a number of clusters of K=4.
}
\usage{data("cora_view1")}
\format{
  The format is:
 num [1:2708, 1:2708] 0 0 0 0 0 0 0 0 1 0 ...
}
\details{
Cora data set includes keyword view, inbound, outbound link view, and citation network view. It takes the form of a sparse matrix. It has 2708 samples and 2708 features.
}
\source{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\references{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\examples{
data(cora_view1); str(cora_view1) 
}
\keyword{datasets}
