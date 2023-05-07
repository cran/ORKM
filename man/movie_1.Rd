\name{movie_1}
\alias{movie_1}
\docType{data}
\title{
The first view of Movie data set.
}
\description{
The first view(keyword view) of Movie data set. The Movie data set contains 2 views, each containing 1878 variables from 617 instances, and the number of clusters to be clustered is K = 17. The number of clusters is large, so it is difficult to cluster. The data set was extracted from IMDb and the main objective was to to find the movie genres, combined from two view matrices.
}
\usage{data("movie_1")}
\format{
  The format is:
 num [1:617, 1:1878] 1 0 0 0 0 0 0 0 0 0 ...
}
\details{
The first view of Movie dataset.
}
\source{
https://lig-membres.imag.fr/grimal/data.html.
}
\references{
C. Grimal. the multi-view movie data set. 2010. URL https://lig-membres.imag.fr/grimal/data.html.
}
\examples{
data(movie_1); str(movie_1) 
}
\keyword{datasets}
