% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/decision2x2.R
\name{serialTables}
\alias{serialTables}
\title{Enumeration of all possible tables}
\usage{
serialTables(N, X0 = rep(0, 4))
}
\arguments{
\item{N}{is integer for number of patients to enter.}

\item{X}{is four values of 2x2 table to start from. Default is c(0,0,0,0)}
}
\value{
Xsm is a list of length N+1, each of which is a list of matrices. The values of matrices are 0 throughout.

Ysm is a list of matrix, each row of which is 4 values of 2x2 table.

M is an utility output.

Z is a list of 4 elements for Aarm&Success, Aarm&Failure, Barm&Succsess, and Barm&Failure; each of which is in the same form with X.
}
\description{
A Function to enumerate all tables for N decisions in a special form
}
\examples{
N <- 20
out <- serialTables(N)
}
\keyword{NA}

