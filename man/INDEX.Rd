\name{INDEX}
\alias{INDEX}
\title{Caculate the indication on the functions
}
\usage{
INDEX(vec1, vec2, method = 0, mybeta = 0)
}
\arguments{
  \item{vec1}{algorithm cluster result}

  \item{vec2}{true cluster result}

  \item{method}{Calculate the selection of indicators.}

  \item{mybeta}{caculate the index }
}
\value{
accuracy
}
\description{
Caculate the indication and the estimator on the functions,method=0 purity;method=1,precision; method=2,recall; method=3, F-score; method=4, RI.
}
\examples{
P1<-c(1,1,1,2,3,2,1);truelabel<-c(1,1,1,2,2,2,3)
INDEX(P1,truelabel,method=0);INDEX(P1,truelabel,method=2)
}
