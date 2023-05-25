\name{cornell_outbound}
\alias{cornell_outbound}
\docType{data}
\title{
The fourth view of   Cornell data set.
}
\description{
 Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data,  Cornell data set,  Texas data set,  Washington data set, and    Wisconsin data set.
}
\usage{data("cornell_outbound")}
\format{
  The format is:
 num [1:195, 1:195] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
 Cornell data set contains four views with a number of clusters of 5. This data set is the fourth view with a sample size of 195 and a number of features of 195.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(cornell_outbound)
## maybe str(cornell_outbound) ; plot(cornell_outbound) ...
}
\keyword{data sets}
