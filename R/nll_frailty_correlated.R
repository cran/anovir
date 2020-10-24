#' Negative log-likelihood function: correlated frailty model
#'
#' Function calculating the negative log-likelihood (nll) for patterns of
#' mortality in infected and uninfected treatments when there is unobserved
#' variation in background mortality and mortality due to infection, where these
#' two sources of variation are positively correlated.
#'
#' This function assumes both sources of unobserved variation follow the Gamma
#' distribution where both have a mean = 1.0, and variances 'theta01' and
#' 'theta02', respectively. The nll function is calculated based on seven
#' parameters; location and scale parameters for background mortality and
#' mortality due to infection, respectively, the unobserved variance of in each
#' source of mortality, and the strength of the positive correlation between
#' them.
#'
#' @param a1,b1 location and scale parameters describing background mortality
#' @param a2,b2 location and scale parameters describing mortality due to
#'   infection
#' @param theta01 parameter describing variance of unobserved variation in
#'   background rate of mortality
#' @param theta02 parameter describing variance of unobserved variation in rate
#'   of mortality due to infection
#' @param rho parameter for strength of positive correlation between theta01 and
#'   theta02; rho > 0
#' @param data name of data frame containing survival data
#' @param time name of data frame column identifying time of event; time > 0
#' @param censor name of data frame column idenifying if event was death (0) or
#'   right-censoring (1)
#' @param infected_treatment name of data frame column identifying if data are
#'   from an infected (1) or uninfected (0) treatment
#' @param d1,d2 names of probability distributions chosen to describe background
#'   mortality and mortality due to infection, respectively; both default to the
#'   Weibull distribution
#' @seealso \code{\link{nll_frailty}}
#' @examples
#'
#' # step #1: parameterise nll function to be passed to 'mle2'
#'     m01_prep_function <- function(
#'       a1 = a1, b1 = b1, a2 = a2, b2 = b2,
#'       theta01 = theta01, theta02 = theta02, rho = rho){
#'         nll_frailty_correlated(
#'           a1 = a1, b1 = b1, a2 = a2, b2 = b2,
#'           theta01 = theta01, theta02 = theta02, rho = rho,
#'           data = data_lorenz,
#'           time = t,
#'           censor = censored,
#'           infected_treatment = g,
#'           d1 = "Gumbel",
#'           d2 = "Gumbel"
#'           )}
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#' #          lower bounds of estimates set to 1e-6
#'   m01 <- mle2(m01_prep_function,
#'     start = list(a1 = 20, b1 = 5, a2 = 20, b2 = 4,
#'                  theta01 = 1, theta02 = 1, rho = 1),
#'     method = "L-BFGS-B",
#'     lower = list(a1 = 1e-6, b1 = 1e-6, a2 = 1e-6, b2 = 1e-6,
#'                  theta01 = 1e-6, theta02 = 1e-6, rho = 1e-6)
#'                  )
#'
#'   summary(m01)
#'
#' # NB no std. errors estimated and estimate of 'theta01' at lower boundary
#'
#' # rerun model with theta01 set at lower boundary
#'     m02 <- mle2(m01_prep_function,
#'       start = list(a1 = 20, b1 = 5, a2 = 20, b2 = 4,
#'                    theta01 = 1, theta02 = 1, rho = 1),
#'       fixed = list(theta01 = 1e-6),
#'       method = "L-BFGS-B",
#'       lower = list(a1 = 1e-6, b1 = 1e-6, a2 = 1e-6,
#'                    b2 = 1e-6, theta02 = 1e-6, rho = 1e-6)
#'       )
#'
#' summary(m02)
#'
#' # NB std. error of 'rho' crosses lower boundary
#'
#' # rerun model with theta01 and rho set at lower limits
#'     m03 <- mle2(m01_prep_function,
#'       start = list(a1 = 20, b1 = 5, a2 = 20, b2 = 4,
#'                    theta01 = 1, theta02 = 1, rho = 1),
#'       fixed = list(theta01 = 1e-6, rho = 1e-6),
#'       method = "L-BFGS-B",
#'       lower = list(a1 = 1e-6, b1 = 1e-6, a2 = 1e-6,
#'                    b2 = 1e-6, theta02 = 1e-6)
#'       )
#'
#'     summary(m03)
#'
#' # result of m03 corresponds with estimates from 'nll_frailty' model,
#'   # i.e., where it is assumed there is no unobserved variation in the
#'   # rate of background mortality and where the gamma distribution describes
#'   # the unobserved variation in virulence
#'
#' # step #1: parameterise nll function to be passed to 'mle2'
#'     m04_prep_function <- function(a1, b1, a2, b2, theta){
#'       nll_frailty(a1, b1, a2, b2, theta,
#'         data = data_lorenz,
#'         time = t,
#'         censor = censored,
#'         infected_treatment = g,
#'         d1 = "Gumbel",
#'         d2 = "Gumbel",
#'         d3 = "Gamma"
#'         )}
#'
#' # step #2: send 'prep_function' to 'mle2' for maximum likelihood estimation
#'     m04 <- mle2(m04_prep_function,
#'       start = list(a1 = 20, b1 = 5, a2 = 20, b2 = 5, theta = 2)
#'       )
#'
#'     summary(m04)
#'
#' # compare estimated values m03 vs. m04
#'     coef(m03) ; coef(m04)
#'
#'
nll_frailty_correlated <- function(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2,
  theta01 = theta01, theta02 = theta02, rho = rho,
  data = data,
  time = time,
  censor = censor,
  infected_treatment = infected_treatment,
  d1 = "Weibull", d2 = "Weibull"
  ){

  pfa1 <- a1
  pfb1 <- b1
  pfa2 <- a2
  pfb2 <- b2

  th1 <- theta01
  th2 <- theta02
  r01 <- rho

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

  # infected population equations

  bl01 <- h1 / (1 + (th1 * H1))
  bl02 <- h2 / (1 + (th2 * H2))

  bl03i   <- (r01 * sqrt(th1) * sqrt(th2))
  bl03ii  <- ((h1 / (1 + th1 * H1)) * H2)
  bl03iii <- ((h2 / (1 + th2 * H2)) * H1)
  bl03iv  <- (1 + th1 * H1 + th2 * H2)

  bl03 <- -(bl03i * ((bl03ii + bl03iii) / bl03iv))

  bl04i   <- (r01 / (sqrt(th1) * sqrt(th2)))
  bl04ii  <- ((1 + th1 * H1)^(-1 / th1))^(-th1)
  bl04iii <- ((1 + th2 * H2)^(-1 / th2))^(-th2)

  bl04 <- -bl04i * log(bl04ii + bl04iii - 1)

  bl05i  <- (1 - (r01 * sqrt(th1) / sqrt(th2)))
  bl05ii <- log((1 + th1 * H1)^(-1 / th1))

  bl05 <- bl05i * bl05ii

  bl06i  <- 1 - (r01 * sqrt(th2) / sqrt(th1))
  bl06ii <- log((1 + th2 * H2)^(-1 / th2))

  bl06 <- bl06i * bl06ii

  # likelihood model

  logl <- 0

  uninfected <- d * log(h1 / (1 + th1 * H1)) + (-1 / th1) * log(1 + th1 * H1)

  infected <- d * log (bl01 + bl02 + bl03) + bl04 + bl05 + bl06

  model <- ifelse(g == 0, uninfected, infected)

  #  print(sum(model))

  if("fq" %in% colnames(data)){
    logl <- sum(data$fq * model)
  } else {
      logl <- sum(model)
  }

  return(-logl)
}

