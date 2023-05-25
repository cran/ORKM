\name{labelcornell}
\alias{labelcornell}
\docType{data}
\title{
True clustering labels for Cornell data set.
}
\description{
Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data, Cornell data set, Texas data set,  Washington data set, and Wisconsin data set.
}
\usage{data("labelcornell")}
\format{
  The format is:
 int [1:195, 1] 1 1 2 3 3 3 2 4 3 3 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "V1"
}
\details{
 Cornell dat aset contains four views with a number of clusters of 5. You can use this true label to calculate your clustering accuracy.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(labelcornell)
## maybe str(labelcornell) ; plot(labelcornell) ...
}
\keyword{data sets}
