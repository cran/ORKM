\name{texas_inbound}
\alias{texas_inbound}
\docType{data}
\title{
The third view of  Texas data set.
}
\description{
 Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data,  Cornell data set,  Texas data set,  Washington data set, and   Wisconsin data set.
}
\usage{data("texas_inbound")}
\format{
  The format is:
 num [1:187, 1:187] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
 Texas data set contains four views with a number of clusters of 5. This data set is the third view with a sample size of 187 and a number of features of 187.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(texas_inbound)
## maybe str(texas_inbound) ; plot(texas_inbound) ...
}
\keyword{data sets}
