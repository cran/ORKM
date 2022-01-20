\name{citesview2}
\alias{citesview2}
\docType{data}
\title{
The second view of Cora dataset.
}
\description{
The second view of Cora dataset.
}
\usage{data("citesview2")}
\format{
  The format is:
Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:49216] 0 0 0 0 0 0 0 0 0 0 ...
  ..@ j       : int [1:49216] 118 125 176 252 351 456 507 521 619 648 ...
  ..@ Dim     : int [1:2] 2708 1433
  ..@ Dimnames:List of 2
  .. ..$ : NULL
  .. ..$ : NULL
  ..@ x       : num [1:49216] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ factors : list()
}
\details{
The second view of Cora dataset.
}
\source{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\references{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\examples{
data(citesview2)
## maybe str(citesview2) ; plot(citesview2) ...
}
\keyword{datasets}
