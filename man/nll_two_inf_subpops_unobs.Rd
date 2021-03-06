% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nll_functions.R
\name{nll_two_inf_subpops_unobs}
\alias{nll_two_inf_subpops_unobs}
\title{Negative log-likelihood function:
two unobserved subpopulations of infected hosts}
\usage{
nll_two_inf_subpops_unobs(
  a1 = a1,
  b1 = b1,
  a2 = a2,
  b2 = b2,
  a3 = a3,
  b3 = b3,
  p1 = p1,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  d1 = "Weibull",
  d2 = "Weibull",
  d3 = "Weibull"
)
}
\arguments{
\item{a1, b1}{location and scale parameters describing background mortality}

\item{a2, b2}{location and scale parameters describing mortality due to
infection in one subpopulation}

\item{a3, b3}{location and scale parameters describing mortality due to
infection in the other subpopulation}

\item{p1}{parameter estimating the proportion of infected hosts in the
first of the two subpopulations; 0 <= p1 <= 1}

\item{data}{name of data frame containing survival data}

\item{time}{name of data frame column identifying time of event; time > 0}

\item{censor}{name of data frame column idenifying if event was
death (0) or right-censoring (1)}

\item{infected_treatment}{name of data frame column identifying if data are
from an infected (1) or uninfected (0) treatment}

\item{d1, d2, d3}{names of probability distributions chosen to describe
background mortality and mortality due to infection in each subpopulation,
respectively; defaults to the Weibull distribution}
}
\description{
Function returning negative log-likelihood (nll) for patterns of mortality in
infected and uninfected treatments when an infected population is assumed
to harbour two distinct subpopulations of hosts experiencing different
virulence.
The nll is based on seven parameters, the location and scale
parameters for background mortality, separate location and scale parameters
for each of the two infected subpopulations, and a parameter estimating the
proportions of the two subpopulations
}
\details{
p1 is the estimated proportion of hosts associated with the location and
scale parameters a2, b2; 0 <= p1 <= 1.

It is assumed the patterns of mortality within each subpopulation
act independently of one another.
}
\examples{
# example using data from Parker et al
    data01 <- data_parker

# step #1: parameterise nll function to be passed to 'mle2'
    m01_prep_function <- function(
      a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3, p1 = p1){
        nll_two_inf_subpops_unobs(
          a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3, p1 = p1,
          data = data01,
          time = t,
          censor = censored,
          infected_treatment = g,
          d1 = "Frechet",
          d2 = "Weibull",
          d3 = "Weibull"
          )}

# step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
    m01 <- mle2(
      m01_prep_function,
      start = list(a1 = 2, b1 = 1,
                   a2 = 2, b2 = 0.3,
                   a3 = 2, b3 = 0.7,
                   p1 = 0.5)
      )

    summary(m01)

}
\seealso{
\code{\link{nll_exposed_infected}}
\code{\link{nll_two_inf_subpops_obs}}
}
