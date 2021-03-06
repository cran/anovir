% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{conf_ints_virulence}
\alias{conf_ints_virulence}
\title{Approximate 95\% confidence intervals for virulence}
\usage{
conf_ints_virulence(
  a2 = a2,
  b2 = b2,
  var_a2 = var_a2,
  var_b2 = var_b2,
  cov_a2b2 = cov_a2b2,
  d2 = "",
  tmax = 21
)
}
\arguments{
\item{a2}{numeric. Estimated value of location parameter describing
mortality due to infection}

\item{b2}{numeric. Estimated value of scale parameter describing
mortality due to infection}

\item{var_a2}{numeric. Estimated variance of location parameter describing
mortality due to infection}

\item{var_b2}{numeric. Estimated variance of scale parameter describing
mortality due to infection}

\item{cov_a2b2}{numeric. Estimated covariance of location and scale
parameters above}

\item{d2}{character. Probability distribution assumed to describe
virulence; Weibull, Gumbel or Fréchet}

\item{tmax}{maximum time virulence will be calculated for. Default value;
tmax = 21}
}
\value{
matrix containing estimates of virulence over time ± approx. 95\%
  confidence intervals
}
\description{
Function calculating the 95\% confidence intervals
for a hazard function based on the variance and covariance
of its location and scale parameters.
}
\details{
The approach is based on the interval being estimated as a complementary
log-log function of the hazard function, h(t), with the variance
of virulence being estimated by the Delta method applied to log(h[t]).
}
\examples{

# the values, variance and covariance of the location and scale parameters
# [a2,a2] describing mortality due to infection were estimated as;
# a2 = 2.5807642
# b2 = 0.1831328
# var_a2 = 0.0008196927
# var_b2 = 0.0010007282
# cov_a2b2 = -0.0003119921

 ci_matrix01 <- conf_ints_virulence(
   a2 = 2.5807642,
   b2 = 0.1831328,
   var_a2 = 0.0008196927,
   var_b2 = 0.0010007282,
   cov_a2b2 = -0.0003119921,
   d2 = "Weibull",
   tmax = 15)

 tail(ci_matrix01)

 plot(ci_matrix01[, 't'], ci_matrix01[, 'h2'],
   type = 'l', col = 'red',
   xlab = 'time', ylab = 'virulence (± 95\% ci)')
   lines(ci_matrix01[, 'lower_ci'], col = 'grey')
   lines(ci_matrix01[, 'upper_ci'], col = 'grey')

}
