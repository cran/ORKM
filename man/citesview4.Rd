\name{citesview4}
\alias{citesview4}
\docType{data}
\title{
The fourth view of Cora dataset.
}
\description{
The fourth view of Cora dataset.
}
\usage{data("citesview4")}
\format{
  The format is:
Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:5429] 402 659 1696 2295 1274 1286 1544 2600 2363 1905 ...
  ..@ j       : int [1:5429] 163 163 163 163 163 163 163 163 163 163 ...
  ..@ Dim     : int [1:2] 2708 2708
  ..@ Dimnames:List of 2
  .. ..$ : NULL
  .. ..$ : NULL
  ..@ x       : num [1:5429] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ factors : list()
}
\details{
The fourth view of Cora dataset.
}
\source{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\references{
http://www.cs.umd.edu/projects/linqs/projects/lbc/
}
\examples{
data(citesview4)
## maybe str(citesview4) ; plot(citesview4) ...
}
\keyword{datasets}
