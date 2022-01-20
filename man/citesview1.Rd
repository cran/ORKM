\name{citesview1}
\alias{citesview1}
\docType{data}
\title{
The first view of Cora dataset.
}
\description{
The first view of Cora dataset.
}
\usage{data("citesview1")}
\format{
  The format is:
Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:10556] 0 0 0 0 0 1 2 2 2 2 ...
  ..@ j       : int [1:10556] 8 14 258 435 544 344 410 471 552 565 ...
  ..@ Dim     : int [1:2] 2708 2708
  ..@ Dimnames:List of 2
  .. ..$ : NULL
  .. ..$ : NULL
  ..@ x       : num [1:10556] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ factors : list()
}
\details{
The first view of Cora dataset.
}
\source{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\references{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\examples{
data(citesview1)
## maybe str(citesview1) ; plot(citesview1) ...
}
\keyword{datasets}
