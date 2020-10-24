# Private function, only for package
# Get z1 as function of distribution type
P_get_z1 <- function(t, a1, b1, d1){
  z1 <- 0
  if(d1 == "Weibull"){
    z1 <- (log(t) - a1) / b1
  } else {
    z1 <- (log(t) - a1) / b1
  }
  return(z1)
}

# Private function, only for package
# Get h1 as function of distribution type
P_get_h1 <- function(t, z1, b1, d1){
  h1 <- 0
  if(d1 == "Weibull"){
    h1 <- (1 / (b1 * t)) * exp(z1)
  } else {
    h1 <- (1 / (b1 * t)) * exp(z1)
  }
  return(h1)
}

# Private function, only for package
# Get S1 as function of distribution type
P_get_S1 <- function(t, z1, d1){
  S1 <- 0
  if(d1 == "Weibull"){
    S1 <- exp(-exp(z1))
  } else {
    S1 <- exp(-exp(z1))
  }
  return(S1)
}


#' Function simulating survival data for nll_basic
#'
#' Function simulating survival data corresponding with assumptions of nll_basic
#'
#' To generate data, the function collects the input variables for the
#' probability distributions chosen to describe the background mortality and
#' mortality due to infection, along with values for their location and scale
#' parameters, respectively.
#' These are used to parameterise cumulative survival functions for the two
#' sources of mortality. These functions are solved for time t, where the
#' value of S(t) is drawn from a random uniform distribution between 0 and 1;
#' this yields values of t corresponding with values of survival
#' drawn at random from the cumulative survival curve.
#'
#' Values of t are rounded up to the nearest integer.
#' This introduces a bias as recorded times of death are later than actual
#' times of death;
#' this bias occurs in real datasets when sampling occurs at discrete intervals.
#' This bias can be partially offset by taking the mid-point of sampling
#' intervals as the time of death.
#' Mid-point times of death are also provided, assuming sampling occurs at
#' each integer interval.
#'
#' In the case of hosts from infected treatments, two potential times of
#' death are calculated;
#' that due to background mortality and that due to infection.
#' The earlier of the two defines the time of death, providing it is not
#' later than the value of t_censor.
#'
#' The value of t_censor defines the time when any individuals remaining alive
#' are right-censored.
#' The act of censoring is assumed to occur after populations have been
#' sampled for mortality, i.e.,
#' if mortality is recorded daily and the last day of an experiment is day 14,
#' the populations are checked for mortality at the beginning of day 14 and
#' any individuals alive after this are classed as censored on day 14.
#' The timing of censoring is 'true' as individuals were known to be alive
#' at the time they were censored. Hence times of censoring do not vary for
#' 'time' or 'mid-time', whereas they do for individuals dying at an unknown
#' time between two consecutive sampling times.
#'
#' @param a1,b1,a2,b2 location and scale parameters describing background
#' mortality and mortality due to infection, respectively; values must be > 0
#' @param n1,n2 size of populations to simulate data; uninfected and infected,
#' respectively; values must be > 0
#' @param t_censor time when any individuals still alive are right-censored;
#' defaults to 1000 if not specified, must be > 0
#' @param d1,d2 probability distribution(s) chosen to describe background
#' mortality and mortality due to infection, respectively;
#' choice of, Weibull, Gumbel, Fréchet
#' @section Warning: The lower tail of Gumbel distribution can yield negative
#' estimates for times of death as the random variable replacing
#' cumulative survival approaches one; S(t) -> 1.
#' @examples
#'   set.seed(34)
#'
#'   sim_data <- sim_data_nll_basic(
#'     a1 = 2.5, b1 = 0.9, n1 = 100, d1 = "Weibull",
#'     a2 = 2.0, b2 = 0.5, n2 = 100, d2 = "Weibull",
#'     t_censor = 14)
#'
#'   head(sim_data, 10)
#'
#'   sim_data$time[sim_data$infected_treatment == 0]
#'
#'   sim_data$time[sim_data$infected_treatment == 1]
#'
#'   # plot histogram for ages at death in infected population
#'       hist(
#'         sim_data$time[(sim_data$infected_treatment == 1 & sim_data$censor == 0)],
#'         breaks = seq(0, 14, 1),
#'         xlim = c(0, 15),
#'         main = 'infected, died',
#'         xlab = 'time of death',
#'         xaxt = 'n')
#'       axis(side = 1, tick = TRUE, at = seq(0, 14, 1))
#'
#'   # estimate location and scale parameters of simulated data with nll_basic
#'       m01_prep_function <- function(a1, b1, a2, b2){
#'         nll_basic(a1, b1, a2, b2,
#'           data = sim_data,
#'           time = time,
#'           censor = censor,
#'           infected_treatment = infected_treatment,
#'           d1 = 'Weibull', d2 = 'Weibull'
#'           )}
#'
#'     m01 <- mle2(m01_prep_function,
#'              start = list(a1 = 2, b1 = 1, a2 = 2, b2 = 1)
#'              )
#'
#'     confint(m01)
#'
#'   # estimate using 'mid_time' instead of 'time'
#'       m02_prep_function <- function(a1, b1, a2, b2){
#'         nll_basic(a1, b1, a2, b2,
#'           data = sim_data,
#'           time = mid_time,
#'           censor = censor,
#'           infected_treatment = infected_treatment,
#'           d1 = 'Weibull', d2 = 'Weibull'
#'           )}
#'
#'       m02 <- mle2(m02_prep_function,
#'              start = list(a1 = 2, b1 = 1, a2 = 2, b2 = 1)
#'              )
#'
#'     confint(m02)
#'
#'     # compare estimates
#'       AICc(m01, m02, nobs = 200)
#'       # for these data, m02 based on 'mid-time' provides a better
#'       # description of the data
sim_data_nll_basic <- function(a1 = a1, b1 = b1, n1 = n1,
                               a2 = a2, b2 = b2, n2 = n2,
                               t_censor = 1000,
                               d1 = "", d2 = ""){


  P_check_inputs_numeric_and_greater_than_zero(
    a1 = a1, b1 = b1, n1 = n1, a2 = a2, b2 = b2, n2 = n2, t_censor = t_censor
  )

  matrix01 <- matrix(0, n1 + n2, 10)

  colnames(matrix01) <- c('infected_treatment', 't1c', 't2c', 't2i',
                          't2ci', 'tci', 'tci_censor',
                          'censor', 'time', 'mid_time')
  i <- 1

  for(i in 1:(n1 + n2)){

    rg01 <- P_generator_random_times(a1, b1, d1)
    rg02 <- P_generator_random_times(a1, b1, d1)
    rg03 <- P_generator_random_times(a2, b2, d2)

    matrix01[i, 1] <- ifelse(i <= n1, 0, 1)
    matrix01[i, 2] <- rg01

    matrix01[i, 3] <- rg02
    matrix01[i, 4] <- rg03
    matrix01[i, 5] <- ifelse(matrix01[i, 3] <= matrix01[i, 4],
                             matrix01[i, 3], matrix01[i, 4])
    matrix01[i, 6] <- ifelse(matrix01[i, 1] == 0,
                             matrix01[i, 2], matrix01[i, 5])
    matrix01[i, 7] <- ifelse(matrix01[i, 6] <= t_censor,
                             matrix01[i, 6], t_censor)
    matrix01[i, 8] <- ifelse(matrix01[i, 6] <= t_censor, 0, 1)
    matrix01[i, 9] <- ceiling(matrix01[i, 7])
    matrix01[i, 10] <- ifelse(matrix01[i, 8] == 1,
                              matrix01[i, 9], matrix01[i, 9] - 0.5)
  }
  output <- as.data.frame(matrix01)

  output[2:7] <- NULL #  remove calculation columns

  return(output)
}



#' Negative log-likelihood function: nll proportional virulence
#'
#' Function assuming hazard functions describing virulence are proportional
#' among infected treatments.
#'
#' The proportional hazards assumption requires, h1(t) / h2(t) = c,
#' where h1(t) and h2(t) are two hazard functions at time t,
#' and c is a constant. In this function the hazard functions for the
#' increased mortality due to infection are assumed to be proportional for
#' different infection treatments within the same experiment.
#'
#' In the default form, one infected treatment is denoted '1' and used as the
#' reference treatment for virulence against which the virulence in a second
#' infected population, '2', is scaled. The default function can be extended
#' to compare multiple treatments against the reference
#' (see vignette: Worked examples II)
#'
#'
#' @param a1,b1 location and scale parameters describing background mortality
#' @param a2,b2 location and scale parameters describing mortality due to
#' infection
#' @param theta a constant describing the proportional relationship
#' of virulence; theta > 0
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column identifying if event was death (0) or
#'   right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#'   from an infected (1) or uninfected (0) treatment
#' @param reference_virulence name of data frame column identifying the
#'   infected treatment to use as a reference when estimating virulence;
#'   value of 1 the reference treatment and 2 for treatment to be compared
#' @param d1,d2 names of probability distributions describing background
#'   mortality and mortality due to infection, respectively; both default to the
#'   Weibull distribution
#' @examples
#' data01 <- data_lorenz
#'
#' # create column 'reference_virulence' with values of 1 and 2 when
#'   # Infectious.dose = 5000 and 160000, respectively, otherwise 0
#'
#' data01$reference_virulence <- ifelse(data01$g == 0, 0,
#'  ifelse(data01$g == 1, ifelse(data01$Infectious.dose == 5000, 1,
#'    ifelse(data01$Infectious.dose == 160000, 2, 0)), 0)
#'  )
#'
#' m01_prep_function <- function(
#'   a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta){
#'     nll_proportional_virulence(
#'       a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
#'       data = data01,
#'       time = t,
#'       censor = censored,
#'       infected_treatment = g,
#'       reference_virulence = reference_virulence,
#'       d1 = 'Gumbel', d2 = 'Weibull')
#'       }
#'
#' m01 <- mle2(m01_prep_function,
#'   start = list(a1 = 23, b1 = 5, a2 = 4, b2 = 0.2, theta = 1)
#'   )
#'
#' summary(m01)
#'
#' # virulence in the high dose treatment is estimated as
#'  # ~ 6x greater than in the low dose treatment
#'

nll_proportional_virulence <- function(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  reference_virulence = reference_virulence,
  d1 = "Weibull", d2 = "Weibull"){

  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2
  pfth <- theta

  t <- data[[deparse(substitute(time))]]
  g <- data[[deparse(substitute(infected_treatment))]]
  d <- data[[deparse(substitute(censor))]] * -1 + 1

  rv <- data[[deparse(substitute(reference_virulence))]]

  z1 <- P_get_zx(t, pfa1, pfb1, d1)
  z2 <- P_get_zx(t, pfa2, pfb2, d2)

  h1 <- P_get_hx(t, z1, pfb1, d1)
  h2 <- P_get_hx(t, z2, pfb2, d2)

  S1 <- P_get_Sx(t, z1, d1)
  S2 <- P_get_Sx(t, z2, d2)

  uninf <- d * log(h1) + log(S1)
  inf_ref <- d * log(h1 + h2) + log(S1) + log(S2)
  inf_else <- d * log(h1 + pfth * h2) + log(S1) + pfth * log(S2)

  logl <- 0

  model <- ifelse(g == 0, uninf,
                  (ifelse(rv == 1, inf_ref,
                          ifelse(rv == 2, inf_else, exp(0))
                  )))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
    logl <- sum(model)
  }
  return(-logl)
}


#' Negative log-likelihood function: basic model on logscale
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
#' on a logscale
#' @param a2,b2 location and scale parameters for mortality due to infection
#' on a logscale
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
#'       nll_basic_logscale(
#'         a1 = a1, b1 = b1, a2 = a2, b2 = b2,
#'         data = data01,
#'         time = t,
#'         censor = censor,
#'         infected_treatment = inf,
#'         d1 = 'Weibull', d2 = 'Weibull')
#'         }
#'
#' # step #2: send 'prep_function' to mle2 for maximum likelihood estimation with
#'         #  starting values specified
#'     m01 <- mle2(m01_prep_function,
#'              start = list(a1 = 1, b1 = 1, a2 = 1, b2 = 1)
#'              )
#'
#'     summary(m01)
#'
nll_basic_logscale <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                      data = data,
                      time = time,
                      censor = censor,
                      infected_treatment = infected_treatment,
                      d1 = "Weibull", d2 = "Weibull"){

  # there are 4 parameter functions; pfa1, pfb1, pfa2, pfb2
  # by default, these functions estimate constants; a1, b1, a2, b2
  # but they can be edited to estimate functions; see examples
  pfa1 <- exp(a1)
  pfb1 <- exp(b1)
  pfa2 <- exp(a2)
  pfb2 <- exp(b2)

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
