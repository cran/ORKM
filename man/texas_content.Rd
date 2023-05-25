\name{texas_content}
\alias{texas_content}
\docType{data}
\title{
The second view of  Texas dataset.
}
\description{
 Webkb data set contains web pages from four universities, with the corresponding clusters categorised as Professor, Student, Program, or Other pages. The data set contains four subsets of data,   Cornell dataset,  Texas dataset,   Washington dataset, and  Wisconsin dataset.
}
\usage{data("texas_content")}
\format{
  The format is:
 num [1:187, 1:1703] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
 Texas data set contains four views with a number of clusters of 5. This data set is the second view with a sample size of 187 and a number of features of 1703.
}
\source{
http://www.cs.cmu.edu/~webkb/
}
\references{
M. Craven, D. DiPasquo, D. Freitag, A. McCallum, T. Mitchell, K. Nigam and S. Slattery. Learning to Extract Symbolic Knowledge from the World Wide Web. Proceedings of the 15th National Conference on Artificial Intelligence (AAAI-98).
}
\examples{
data(texas_content)
## maybe str(texas_content) ; plot(texas_content) ...
}
\keyword{data sets}
