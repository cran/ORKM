\name{citesview3}
\alias{citesview3}
\docType{data}
\title{
The third view of Cora dataset.
}
\description{
The third view of Cora dataset.
}
\usage{data("citesview3")}
\format{
  The format is:
Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:5429] 163 163 163 163 163 163 163 163 163 163 ...
  ..@ j       : int [1:5429] 402 659 1696 2295 1274 1286 1544 2600 2363 1905 ...
  ..@ Dim     : int [1:2] 2708 2708
  ..@ Dimnames:List of 2
  .. ..$ : NULL
  .. ..$ : NULL
  ..@ x       : num [1:5429] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ factors : list()
}
\details{
The third view of Cora dataset.
}
\source{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\references{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\examples{
data(citesview3)
## maybe str(citesview3) ; plot(citesview3) ...
}
\keyword{datasets}
