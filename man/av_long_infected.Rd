% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{av_long_infected}
\alias{av_long_infected}
\title{Average longevity: estimate for infected hosts}
\usage{
av_long_infected(a1, b1, a2, b2, d1 = "", d2 = "")
}
\arguments{
\item{a1, b1}{numeric: location & scale parameters for background mortality,
respectively}

\item{a2, b2}{numeric: location & scale parameters for mortality due to
infection, respectively}

\item{d1, d2}{character: probability distributions to describe background
mortality and mortality due to infection, respectively}
}
\value{
a vector
}
\description{
Calculates expected longevity of infected hosts due to background mortality
and mortality due to infection
}
\details{
The expected average longevity is calculated as the integral
from zero to infinity for the product of the cumulative survival
functions for background mortality and mortality due to infection,
given values of a1, b1, d1, a2, b2, d2
}
\examples{
  av_long_infected(
    a1 = 3.0, b1 = 0.6, d1 = "Weibull",
    a2 = 2.5, b2 = 0.5, d2 = "Frechet") # 12.156
}
\seealso{
\code{\link{av_long_uninfected}}
}
