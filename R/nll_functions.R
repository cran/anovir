# NB replaced '≤' with '<=' as unicode '\u2264' does not work

#' Negative log-likelihood function: basic model
#'
#' Function returning the negative log-likelihood (nll) for the 'basic' relative
#' survival model, given the functions' parameters and the observed data.
#'
#' By deafult, this function takes arguments for location and scale parameters
#' named; a1, b1, a2, b2. These parameters are components of survival functions
#' describing patterns of background mortality and mortality due to infection.
#' The particular form of these survival functions depends on the probability
#' distributions chosen to describe each source of mortality; d1, d2. The
#' function also takes arguments directing it to the data to be analysed and how
#' they are labelled; data, time, censor, etc.
#'
#' The nll returned by the function depends on the numerical values of the
#' location and scale parameters, which determine how well the likelihood model
#' describes the observed data. Maximum likelihood estimation functions, e.g.,
#' \code{mle2} of the package \code{bbmle}, can be used find values of the
#' location and scale parameters minimising the model's nll. The resulting
#' maximum likelihood estimates can be used to describe host mortality due to
#' background mortality and mortality due to infection, including the pathogen's
#' virulence.
#'
#' The model assumes all the individuals in the infected population are
#' infected. It is also assumes infections are homogeneous, i.e., each
#' infection has the same influence on host survival. Consequently a single
#' hazard function, with a single pair of values for its location and scale
#' parameters, can be used to describe the pattern of mortality due to infection
#' for the infected population as a whole.
#'
#'
#' @param a1,b1 location and scale parameters for background mortality
#' @param a2,b2 location and scale parameters for mortality due to infection
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column identifying if event was death (0) or
#'   right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#'   from an infected (1) or uninfected (0) treatment
#' @param d1,d2 names of probability distributions describing background
#'   mortality and mortality due to infection, respectively; both default to the
#'   Weibull distribution
#' @examples
#' # prepare subset of 'data_blanford'; treatments 'cont' and 'Bb06' of Block 3
#'
#'   data01 <- subset(data_blanford,
#'     (data_blanford$block == 3) & (
#'       (data_blanford$treatment == 'cont') |
#'          (data_blanford$treatment == 'Bb06')) &
#'            (data_blanford$day > 0))
#'
#'  head(data01, 4)
#'
#' # step #1: 'prep function' linking 'nll_basic' to data
#'   # and identifying parameters to estimate
#'     m01_prep_function <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2){
#'       nll_basic(
#'         a1 = a1, b1 = b1, a2 = a2, b2 = b2,
#'         data = data01,
#'         time = t,
#'         censor = censor,
#'         infected_treatment = inf,
#'         d1 = 'Weibull', d2 = 'Weibull')
#'         }
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#'   #  starting values specified as
#'     m01 <- mle2(m01_prep_function,
#'              start = list(a1 = 2, b1 = 0.5, a2 = 2.5, b2 = 0.25)
#'              )
#'
#'     summary(m01)
#'
nll_basic <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                      data = data,
                      time = time,
                      censor = censor,
                      infected_treatment = infected_treatment,
                      d1 = "Weibull", d2 = "Weibull"){

  # there are 4 parameter functions; pfa1, pfb1, pfa2, pfb2
  # by default, these functions estimate constants; a1, b1, a2, b2
  # but they can be edited to estimate functions; see examples
  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)

  logl <- 0

  model <- (d * log(h1 + g * h2) + log(S1) + g * log(S2))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
    logl <- sum(model)
  }
  return(-logl)
}



#' Negative log-likelihood function: control data only
#'
#' Function returning negative log-likelihood (nll) for data in a control
#' treatment.
#'
#' This function returns the nll based on two parameters, the location and scale parameters
#' used to describe background mortality.
#'
#' @param a1,b1 location and scale parameters for background mortality
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0) or
#'   right-censoring (1)
#' @param d1 name of probability distribution describing background
#'   mortality. Choice of; 'Weibull', 'Gumbel', 'Fréchet'; defaults to the
#'   Weibull distribution
#' @return numeric
#' @seealso \code{\link{nll_basic}}
#' @examples
#' # prepare a subset of the Blanford data for analysis
#'   data01 <- subset(data_blanford,
#'     (data_blanford$block == 3) &
#'       (data_blanford$treatment == 'cont') &
#'         (data_blanford$day > 0))
#'
#' # check data frame for names of columns
#'     head(data01)
#'
#' # step #1: 'prep function' linking 'nll_controls' to data
#'   # and identifying parameters to estimate
#'     m01_prep_function <- function(a1 = a1, b1 = b1){
#'       nll_controls(
#'         a1 = a1, b1 = b1,
#'         data = data01,
#'         time = t,
#'         censor = censor,
#'         d1 = 'Weibull'
#'         )}
#'
#' # step #2: send 'prep_function' to mle2 for maximum likelihood estimation
#'   # specifying starting values
#'     m01 <- mle2(m01_prep_function,
#'              start = list(a1 = 2, b1 = 0.5)
#'              )
#'
#'     summary(m01)
#'
#'
nll_controls <- function(
  a1 = a1,
  b1 = b1,
  data = data,
  time = time,
  censor = censor,
  d1 = "Weibull"){

  pfa1 <- a1
  pfb1 <- b1

  t <- data[[deparse(substitute(time))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  h1 <- P_get_hx(t, z1, pfb1, d1)
  S1 <- P_get_Sx(t, z1, d1)

  logl <- 0

  model <- (d * log(h1) + log(S1))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
    logl <- sum(model)
  }
  return(-logl)
}


#' Negative log-likelihood function: exposed-infected
#'
#' Function returning negative log-likelihood (nll) for patterns of mortality in
#' infected and control treatments, where the infected population harbours an
#' unobserved proportion of hosts that were exposed to infection, but did not
#' become infected.
#'
#' This function returns the nll based on five parameters, the location and
#' scale parameters for background mortality and mortality due to infection,
#' respectively, plus a parameter for the proportion of hosts that became
#' infected when exposed to infection.
#'
#' @param a1,b1 location and scale parameters for background mortality
#' @param a2,b2 location and scale parameters for mortality due to infection
#' @param p1 unobserved proportion of hosts exposed to infection and infected; 0
#'   <= p1 <= 1
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0) or
#'   right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#'   from an infected (1) or uninfected (0) treatment
#' @param d1,d2 names of probability distributions chosen to describe background
#'   mortality and mortality due to infection, respectively'; both default to
#'   the Weibull distribution
#' @seealso \code{\link{nll_two_inf_subpops_obs}}
#'   \code{\link{nll_two_inf_subpops_unobs}}
#' @return numeric
#' @examples
#' # check column names in head of data frame with data to analyse
#'     head(data_parker)
#'
#' # step #1: prepare nll function for analysis
#'     m01_prep_function <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2, p1 = p1){
#'       nll_exposed_infected(
#'         a1 = a1, b1 = b1, a2 = a2, b2 = b2, p1 = p1,
#'         data = data_parker,
#'         time = t,
#'         censor = censored,
#'         infected_treatment = g,
#'         d1 = "Frechet",
#'         d2 = "Weibull")
#'         }
#'
#' # step #2: send 'prep_function' to mle2 for maximum likelihood estimation
#'     m01 <- mle2(m01_prep_function,
#'              start = list(a1 = 2.5, b1 = 1, a2 = 2, b2 = 0.5, p1 = 0.5)
#'              )
#'
#'     summary(m01)
#'
#' # model setting lower & upper bounds to parameter estimates
#'   # including 0 < p1 < 1
#'     m02 <- mle2(m01_prep_function,
#'              start = list(a1 = 2.5, b1 = 1.2, a2 = 1.9, b2 = 0.16, p1 = 0.48),
#'              method = "L-BFGS-B",
#'              lower = c(a1 = 0, b1 = 0, a2 = 0, b2 = 0, p1 = 0),
#'              upper = c(a1 = Inf, b1 = Inf, a2 = Inf, b2 = Inf, p1 = 1),
#'              )
#'
#'     summary(m02)
#'
#'
nll_exposed_infected <- function(
  a1 = a1,
  b1 = b1,
  a2 = a2,
  b2 = b2,
  p1 = p1,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  d1 = "Weibull",
  d2 = "Weibull"
){

  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2
  pfp1 <- p1

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * - 1 + 1

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)

  f1 <- P_get_fx(t, z1, pfb1, d1)
  f2 <- P_get_fx(t, z2, pfb2, d2)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)

  p <- pfp1

  Sobsinf <- p * (S1 * S2) + (1 - p) * S1
  fobsinf <- p * (f1 * S2 + f2 * S1) + (1 - p) * f1
  hobsinf <- fobsinf / Sobsinf

  uninf.treat <- d * log(h1) + log(S1)
  inf.treat <- d * log(hobsinf) + log(Sobsinf)

  logl <- 0

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * (ifelse(g == 0, uninf.treat, inf.treat)))
  } else {
    logl <- sum(ifelse(g == 0, uninf.treat, inf.treat))
  }
  return(-logl)
}


#' Negative log-likelihood function: frailty
#'
#' Function calculating negative log-likelihood (nll) for the observed patterns
#' of mortality in infected and uninfected treatments when it assumed there is
#' unobserved variation in virulence.
#'
#' The function assumes the unobserved variation in the rate of mortality due to
#' infection is continuously distributed and follows either the gamma
#' distribution or the inverse Gaussian distribution, with mean = 1 and variance
#' = theta.
#'
#' The nll is based on five parameter functions for the location and
#' scale parameters for background mortality and mortality due to infection,
#' respectively, plus a parameter estimating the variance of the unobserved
#' variation in virulence.
#'
#' @param a1,b1 location and scale parameters describing background mortality
#' @param a2,b2 location and scale parameters describing mortality due to
#'   infection
#' @param theta parameter describing variance of unobserved variation in
#'   virulence
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0) or
#'   right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#'   from an infected (1) or uninfected (0) treatment
#' @param d1,d2 names of probability distributions chosen to describe background
#'   mortality and mortality due to infection, respectively; both default to the
#'   Weibull distribution
#' @param d3 name of probability distribution chosen to describe unobserved
#'   frailty; choice of 'gamma' or 'inverse Gaussian'
#' @return numeric
#' @examples
#' ### Example 1: unobserved variation in virulence described by gamma distribution
#'
#' # step #1: parameterise nll function to be passed to 'mle2'
#'     m01_prep_function <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta){
#'       nll_frailty(
#'         a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
#'         data = data_lorenz,
#'         time = t,
#'         censor = censored,
#'         infected_treatment = g,
#'         d1 = "Gumbel",
#'         d2 = "Weibull",
#'         d3 = "Gamma"
#'         )}
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#'     m01 <- mle2(
#'       m01_prep_function,
#'       start = list(a1 = 20, b1 = 5, a2 = 3, b2 = 0.1, theta = 2)
#'       )
#'
#'     summary(m01)
#'
#' ### Example 2: unobserved variation in virulence described by inverse Gaussian distribution
#'
#'     m02_prep_function <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta){
#'       nll_frailty(
#'         a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
#'         data = data_lorenz,
#'         time = t,
#'         censor = censored,
#'         infected_treatment = g,
#'         d1 = "Gumbel",
#'         d2 = "Weibull",
#'         d3 = "Inverse Gaussian"
#'         )}
#'
#'     m02 <- mle2(
#'       m02_prep_function,
#'         start = list(a1 = 20, b1 = 5, a2 = 3, b2 = 0.1, theta = 2)
#'         )
#'
#'     summary(m02)
#'
#' # compare model estimates by AICc
#'     AICc(m01, m02, nobs = 256)
#'
nll_frailty <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
                        data = data,
                        time = time,
                        censor = censor,
                        infected_treatment = infected_treatment,
                        d1 = "Weibull",
                        d2 = "Weibull",
                        d3 = ""
){

  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2

  pftheta <- theta

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)

  H1 <- P_get_Hx(t, z1, d1)
  H2 <- P_get_Hx(t, z2, d2)

  theta <- pftheta

  logl <- 0

  if((d3 == "Gamma") | (d3 == "gamma")){
    model <- (d * log(h1 + g * (h2 / (1 + theta * H2))) + log(S1) + g * -(1 / theta) * log(1 + theta * H2))
  } else if ((d3 == "inverse Gaussian") | (d3 == "inverse gaussian") | (d3 == "Inverse Gaussian")) {
    model <- (d * log(h1 + g * (h2 / (1 + 2 * theta * H2)^0.5)) + log(S1) + g * (1 / theta) * (1 - (1 + 2 * theta * H2)^0.5))
  } else {
    model <- 1
  }

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
    logl <- sum(model)
  }
  return(-logl)
}



#' Negative log-likelihood function: frailty variables on logscale
#'
#' This negative log-likelihood (nll) function is the same as 'nll_frailty', except it assumes the variables to estimate are on a logscale.

#' @param log.a1,log.b1 location and scale parameters for background mortality,
#' on a logscale
#' @param log.a2,log.b2 location and scale parameters for mortality
#' due to infection, on a logscale
#' @param log.theta parameter describing variance of the unobserved variation
#' in virulence, on a logscale
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0)
#' or right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data
#' are from an infected (1) or uninfected (0) treatment
#' @param d1,d2 names of probability distributions chosen to describe
#' background mortality and mortality due to infection, respectively;
#' both default to the Weibull distribution
#' @param d3 name of probability distribution chosen to describe
#' unobserved frailty; choice of 'gamma' or 'inverse Gaussian'
#' @return numeric
#' @seealso \code{\link{nll_frailty}}
#' @examples
#'
#' ### Example 1: unobserved variation in virulence with gamma distribution
#'
#' # step #1: parameterise nll function to be passed to 'mle2'
#'     m01_prep_function <- function(
#'       log.a1 = log.a1, log.b1 = log.b1,
#'       log.a2 = log.a2, log.b2 = log.b2,
#'       log.theta = log.theta){
#'         nll_frailty_logscale(
#'           log.a1 = log.a1, log.b1 = log.b1,
#'           log.a2 = log.a2, log.b2 = log.b2,
#'           log.theta = log.theta,
#'           data = data_lorenz,
#'           time = t,
#'           censor = censored,
#'           infected_treatment = g,
#'           d1 = "Gumbel", d2 = "Weibull", d3 = "Gamma"
#'           )}
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#'     m01 <- mle2(
#'       m01_prep_function,
#'       start = list(
#'         log.a1 = 3, log.b1 = 1.5, log.a2 = 0.7, log.b2 = -0.7, log.theta = 1
#'         )
#'       )
#'
#'   summary(m01)
#'
#'   exp(coef(m01))
#'
nll_frailty_logscale <- function(log.a1 = log.a1, log.b1 = log.b1,
                                 log.a2 = log.a2, log.b2 = log.b2,
                                 log.theta = log.theta,
                                 data = data,
                                 time = time,
                                 censor = censor,
                                 infected_treatment = infected_treatment,
                                 d1 = "Weibull",
                                 d2 = "Weibull",
                                 d3 = ""
){

  pflog.a1 <- log.a1
  pflog.b1 <- log.b1
  pflog.a2 <- log.a2
  pflog.b2 <- log.b2

  log.theta <- log.theta

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  z1 <- P_get_zx_logscale(t, pflog.a1, pflog.b1, d1)
  z2 <- P_get_zx_logscale(t, pflog.a2, pflog.b2, d2)

  h1 <- P_get_hx_logscale(t, z1, pflog.b1, d1)
  h2 <- P_get_hx_logscale(t, z2, pflog.b2, d2)

  S1 <- P_get_Sx_logscale(t, z1, d1)
  S2 <- P_get_Sx_logscale(t, z2, d2)

  H1 <- P_get_Hx(t, z1, d1)
  H2 <- P_get_Hx(t, z2, d2)

  logl <- 0

  if((d3 == "Gamma") | (d3 == "gamma")){
    model <- (d * log(h1 + g * (h2 / (1 + exp(log.theta) * H2))) + log(S1) + g * -(1 / exp(log.theta)) * log(1 + exp(log.theta) * H2))
  } else if ((d3 == "Inverse Gaussian") | (d3 == "inverse gaussian")) {
    model <- (d * log(h1 + g * (h2 / (1 + 2 * exp(log.theta) * H2)^0.5)) + log(S1) + g * (1 / exp(log.theta)) * (1 - (1 + 2 * exp(log.theta) * H2)^0.5))
  } else {
    model <- 1
  }

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
    logl <- sum(model)
  }
  return(-logl)
}




#' Negative log-likelihood function: frailty shared
#'
#' Function calculating negative log-likelihood (nll) for patterns of mortality
#' in infected and uninfected treatments where unobserved variation is assumed
#' to act equally on background mortality and mortality due to infection.
#'
#' This function assumes unobserved variation acting on both the background rate
#' of mortality and the rate of mortality due to infection is continuously
#' distributed and follows the gamma distribution, with mean = 1.0 and variance
#' = theta. The function returns the nll based on five parameters; the location
#' and scale parameters for background mortality and mortality due to infection,
#' plus the parameter describing the variance of the unobserved variation.
#'
#' @param a1,b1 location and scale parameters for background mortality
#' @param a2,b2 location and scale parameters for mortality due to infection
#' @param theta parameter describing variance of unobserved variation acting on
#'   mortality rates
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0) or
#'   right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#'   from an infected (1) or uninfected (0) treatment
#' @param d1,d2 names of probability distributions chosen to describe background
#'   mortality and mortality due to infection, respectively; both default to the
#'   Weibull distribution
#' @return numeric
#' @seealso \code{\link{nll_frailty}}
#' @examples
#'
#' # step #1: prepare nll function for analysis
#'   m01_prep_function <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta){
#'     nll_frailty_shared(a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
#'       data = data_lorenz,
#'       time = t,
#'       censor = censored,
#'       infected_treatment = g,
#'       d1 = "Gumbel", d2 = "Gumbel"
#'       )}
#'
#' # step #2: send 'prep_function' to mle2 for maximum likelihood estimation,
#'   # specifying starting values
#'   m01 <- mle2(m01_prep_function,
#'             start = list(a1 = 23, b1 = 5, a2 = 10, b2 = 1, theta = 1),
#'             method = "Nelder-Mead",
#'             control = list(maxit = 5000)
#'             )
#'
#'   summary(m01)
#'
nll_frailty_shared <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
                               data = data,
                               time = time,
                               censor = censor,
                               infected_treatment = infected_treatment,
                               d1 = "", d2 = ""){

  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)

  H1 <- P_get_Hx(t, z1, d1)
  H2 <- P_get_Hx(t, z2, d2)

  logl <- 0

  model <- d * log((h1 + g * h2) / (1 + theta * (H1 + g * H2))) + log((1 + theta * (H1 + g * H2))^(-1 / theta))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
    logl <- sum(model)
  }
  return(-logl)
}

#' Negative log-likelihood function:
#' two observed subpopulations of infected hosts
#'
#' Function returning negative log-likelihood (nll) for patterns of mortality in
#' infected and uninfected treatments when an infected population harbours two
#' identified, or 'observed', subpopulations of hosts experiencing
#' different patterns of virulence,
#' e.g. with/without visible signs of infection.
#'
#' The nll is based on six parameters, the location and scale
#' parameters for background mortality,
#' plus separate location and scale parameters for each of
#' the two infected subpopulations.
#'
#' It is assumed the patterns of mortality within each subpopulation
#' act independently of one another.
#'
#' @param a1,b1 location and scale parameters describing background mortality
#' @param a2,b2 location and scale parameters describing mortality due to
#' infection in one subpopulation
#' @param a3,b3 location and scale parameters describing mortality due to
#' infection in the other subpopulation
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0)
#' or right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data
#' are from an infected (1) or uninfected (0) treatment
#' @param d1,d2,d3 names of probability distributions chosen to describe
#' background mortality and mortality due to infection, respectively;
#' each defaults to the Weibull distribution
#' @param infsubpop name of data frame column identifying
#' the two subpopulations of infected hosts; '1' or '2'
#' @seealso \code{\link{nll_exposed_infected}}
#' \code{\link{nll_two_inf_subpops_unobs}}
#' @examples
#' # example using data from Parker et al
#'     data01 <- data_parker
#'
#' # create column 'infsubpop' in data01, fill with '0'
#'     data01$infsubpop <- 0
#'
#' # infsubpop = '1' for individuals in infected treatments (g == 1)
#'   # with visible signs of sporulation (Sporulation = 1)
#' # infsubpop = '2' for individuals in infected treatments (g == 1)
#'   # with no visible signs of sporulation (Sporulation = 0)
#'     data01$infsubpop[data01$g == 1 & data01$Sporulation == 1] <- 1
#'     data01$infsubpop[data01$g == 1 & data01$Sporulation == 0] <- 2
#'
#'     head(data01)
#'
#' # step #1: parameterise nll function to be passed to 'mle2'
#'     m01_prep_function <- function(
#'       a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3){
#'         nll_two_inf_subpops_obs(
#'           a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3,
#'           data = data01,
#'           time = t,
#'           censor = censored,
#'           infected_treatment = g,
#'           d1 = "Frechet",
#'           d2 = "Weibull",
#'           d3 = "Weibull",
#'           infsubpop = infsubpop
#'         )}
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#'     m01 <- mle2(
#'       m01_prep_function,
#'       start = list(a1 = 3, b1 = 1, a2 = 2, b2 = 0.5, a3 = 2, b3 = 0.5)
#'       )
#'
#'     summary(m01)
#'
nll_two_inf_subpops_obs <- function(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  d1 = "Weibull",
  d2 = "Weibull",
  d3 = "Weibull",
  infsubpop = infsubpop){
  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2
  pfa3 <- a3
  pfb3 <- b3

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  infsubpop <- data[[deparse(substitute(infsubpop))]]

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)
  z3 <- P_get_zx(t, pfa3, pfb3, d3)

  f1 <- P_get_fx(t, z1, pfb1, d1)
  f2 <- P_get_fx(t, z2, pfb2, d2)
  f3 <- P_get_fx(t, z3, pfb3, d3)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)
  h3 <- P_get_hx(t, z3, pfb3, d3)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)
  S3 <- P_get_Sx(t, z3, d3)

  S1S2 <- S1 * S2
  S1S3 <- S1 * S3

  fobs1 <- (f1 * S2) + (f2 * S1)
  fobs2 <- (f1 * S3) + (f3 * S1)

  hobs1 <- fobs1 / S1S2
  hobs2 <- fobs2 / S1S3

  uninfected_hosts <- d * log(h1) + log(S1)
  infsubpop1 <- d * log(h1 + h2) + log(S1) + log(S2)
  infsubpop2 <- d * log(h1 + h3) + log(S1) + log(S3)

  logl <- 0

  model <- (ifelse(g == 0, uninfected_hosts,
                   ifelse(infsubpop == 1, infsubpop1, infsubpop2)))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
      logl <- sum(model)
  }

  return(-logl)
}



#' Negative log-likelihood function:
#' two unobserved subpopulations of infected hosts
#'
#' Function returning negative log-likelihood (nll) for patterns of mortality in
#' infected and uninfected treatments when an infected population is assumed
#' to harbour two distinct subpopulations of hosts experiencing different
#' virulence.
#' The nll is based on seven parameters, the location and scale
#' parameters for background mortality, separate location and scale parameters
#' for each of the two infected subpopulations, and a parameter estimating the
#' proportions of the two subpopulations
#'
#' p1 is the estimated proportion of hosts associated with the location and
#' scale parameters a2, b2; 0 <= p1 <= 1.
#'
#' It is assumed the patterns of mortality within each subpopulation
#' act independently of one another.
#' @param a1,b1 location and scale parameters describing background mortality
#' @param a2,b2 location and scale parameters describing mortality due to
#' infection in one subpopulation
#' @param a3,b3 location and scale parameters describing mortality due to
#' infection in the other subpopulation
#' @param p1 parameter estimating the proportion of infected hosts in the
#' first of the two subpopulations; 0 <= p1 <= 1
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was
#' death (0) or right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#' from an infected (1) or uninfected (0) treatment
#' @param d1,d2,d3 names of probability distributions chosen to describe
#' background mortality and mortality due to infection in each subpopulation,
#' respectively; defaults to the Weibull distribution
#' @seealso \code{\link{nll_exposed_infected}}
#' \code{\link{nll_two_inf_subpops_obs}}
#' @examples
#' # example using data from Parker et al
#'     data01 <- data_parker
#'
#' # step #1: parameterise nll function to be passed to 'mle2'
#'     m01_prep_function <- function(
#'       a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3, p1 = p1){
#'         nll_two_inf_subpops_unobs(
#'           a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3, p1 = p1,
#'           data = data01,
#'           time = t,
#'           censor = censored,
#'           infected_treatment = g,
#'           d1 = "Frechet",
#'           d2 = "Weibull",
#'           d3 = "Weibull"
#'           )}
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#'     m01 <- mle2(
#'       m01_prep_function,
#'       start = list(a1 = 2, b1 = 1,
#'                    a2 = 2, b2 = 0.3,
#'                    a3 = 2, b3 = 0.7,
#'                    p1 = 0.5)
#'       )
#'
#'     summary(m01)
#'
nll_two_inf_subpops_unobs <- function(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2, a3 = a3, b3 = b3, p1 = p1,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  d1 = "Weibull",
  d2 = "Weibull",
  d3 = "Weibull"
  ){

  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2
  pfa3 <- a3
  pfb3 <- b3

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)
  z3 <- P_get_zx(t, pfa3, pfb3, d3)

  f1 <- P_get_fx(t, z1, pfb1, d1)
  f2 <- P_get_fx(t, z2, pfb2, d2)
  f3 <- P_get_fx(t, z3, pfb3, d3)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)
  h3 <- P_get_hx(t, z3, pfb3, d3)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)
  S3 <- P_get_Sx(t, z3, d3)

  p <- p1

  Sobsinf <- p * (S1 * S2) + (1 - p) * (S1 * S3)
  fobsinf <- p * (f1 * S2 + f2 * S1) + (1 - p) * (f1 * S3 + f3 * S1)
  hobsinf <- fobsinf / Sobsinf

  uninf.treat <- d * log(h1) + log(S1)
  inf.treat <- d * log(hobsinf) + log(Sobsinf)

  logl <- 0

  model <- (ifelse(g == 0, uninf.treat, inf.treat))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model )
  } else {
    logl <- sum(model)
  }
  return(-logl)
}


