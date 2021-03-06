% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nll_functions.R
\name{nll_two_inf_subpops_obs}
\alias{nll_two_inf_subpops_obs}
\title{Negative log-likelihood function:
two observed subpopulations of infected hosts}
\usage{
nll_two_inf_subpops_obs(
  a1 = a1,
  b1 = b1,
  a2 = a2,
  b2 = b2,
  a3 = a3,
  b3 = b3,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  d1 = "Weibull",
  d2 = "Weibull",
  d3 = "Weibull",
  infsubpop = infsubpop
)
}
\arguments{
\item{a1, b1}{location and scale parameters describing background mortality}

\item{a2, b2}{location and scale parameters describing mortality due to
infection in one subpopulation}

\item{a3, b3}{location and scale parameters describing mortality due to
infection in the other subpopulation}

\item{data}{name of data frame containing survival data}

\item{time}{name of data frame column identifying time of event; time > 0}

\item{censor}{name of data frame column idenifying if event was death (0)
or right-censoring (1)}

\item{infected_treatment}{name of data frame column identifying if data
are from an infected (1) or uninfected (0) treatment}

\item{d1, d2, d3}{names of probability distributions chosen to describe
background mortality and mortality due to infection, respectively;
each defaults to the Weibull distribution}

\item{infsubpop}{name of data frame column identifying
the two subpopulations of infected hosts; '1' or '2'}
}
\description{
Function returning negative log-likelihood (nll) for patterns of mortality in
infected and uninfected treatments when an infected population harbours two
identified, or 'observed', subpopulations of hosts experiencing
different patterns of virulence,
e.g. with/without visible signs of infection.
}
\details{
The nll is based on six parameters, the location and scale
parameters for background mortality,
plus separate location and scale parameters for each of
the two infected subpopulations.

It is assumed the patterns of mortality within each subpopulation
act independently of one another.
}
\examples{
# example using data from Parker et al
    data01 <- data_parker

# create column 'infsubpop' in data01, fill with '0'
    data01$infsubpop <- 0

# infsubpop = '1' for individuals in infected treatments (g == 1)
  # with visible signs of sporulation (Sporulation = 1)
# infsubpop = '2' for individuals in infected treatments (g == 1)
  # with no visible signs of sporulation (Sporulation = 0)
    data01$infsubpop[data01$g == 1 & data01$Sporulation == 1] <- 1
    data01$infsubpop[data01$g == 1 & data01$Sporulation == 0] <- 2

    head(data01)

# step #1: parameterise nll function to be passed to 'mle2'
    m01_prep_function <- function(
      a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3){
        nll_two_inf_subpops_obs(
          a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3,
          data = data01,
          time = t,
          censor = censored,
          infected_treatment = g,
          d1 = "Frechet",
          d2 = "Weibull",
          d3 = "Weibull",
          infsubpop = infsubpop
        )}

# step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
    m01 <- mle2(
      m01_prep_function,
      start = list(a1 = 3, b1 = 1, a2 = 2, b2 = 0.5, a3 = 2, b3 = 0.5)
      )

    summary(m01)

}
\seealso{
\code{\link{nll_exposed_infected}}
\code{\link{nll_two_inf_subpops_unobs}}
}
