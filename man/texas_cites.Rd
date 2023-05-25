\name{texas_cites}
\alias{texas_cites}
\docType{data}
\title{
The first view of  Texas data set.
}
\description{
 Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data,  Cornell data set, Texas data set,  Washington data set, and Wisconsin data set.
}
\usage{data("texas_cites")}
\format{
  The format is:
 num [1:187, 1:187] 0 1 1 1 0 1 1 0 1 0 ...
}
\details{
 Texas data set contains four views with a number of clusters of 5. This data set is the first view with a sample size of 187 and a number of features of 187.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(texas_cites)
## maybe str(texas_cites) ; plot(texas_cites) ...
}
\keyword{data sets}
