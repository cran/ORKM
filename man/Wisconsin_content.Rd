\name{Wisconsin_content}
\alias{Wisconsin_content}
\docType{data}
\title{
The second view of Wisconsin data set.
}
\description{
 Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data, Cornell data set,  Texas data set,  Washington data set, and  Wisconsin data set.
}
\usage{data("Wisconsin_content")}
\format{
  The format is:
 num [1:265, 1:1703] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
 Wisconsin data set contains four views with a number of clusters of 5. This data set is the second view with a sample size of 265 and a number of features of 1703.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(Wisconsin_content)
## maybe str(Wisconsin_content) ; plot(Wisconsin_content) ...
}
\keyword{data sets}
