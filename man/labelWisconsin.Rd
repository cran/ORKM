\name{labelWisconsin}
\alias{labelWisconsin}
\docType{data}
\title{
True clustering labels for  Wisconsin data set.
}
\description{
  Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data,  Cornell dataset,  Texas dataset,   Washington dataset, and   Wisconsin data set.
}
\usage{data("labelWisconsin")}
\format{
  The format is:
 num [1:265] 1 2 3 3 1 1 1 1 1 1 ...
}
\details{
 Wisconsin data set contains four views with a number of clusters of 5. You can use this true label to calculate your clustering accuracy.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(labelWisconsin)
## maybe str(labelWisconsin) ; plot(labelWisconsin) ...
}
\keyword{data sets}
