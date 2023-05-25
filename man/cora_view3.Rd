\name{cora_view3}
\alias{cora_view3}
\docType{data}
\title{
The third view of Cora data set.
}
\description{
This data matrix is the third view of Cora data set. It called the inbound link view and the form of a sparse matrix. It has 2708 samples and 2708 features.  Cora data set is a multi-view data set of machine learning papers with 4 views, a sample size of nearly 3000 and a number of features of 1500, with a number of clusters of K=4.
}
\usage{data("cora_view3")}
\format{
  The format is:
 num [1:2708, 1:2708] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
The third view of Cora data set.
}
\source{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\references{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\examples{
data(cora_view3); str(cora_view3) 
}
\keyword{data sets}
