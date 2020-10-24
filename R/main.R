


#' Average longevity: estimate for uninfected hosts
#'
#' Calculates expected average longevity due only to background mortality
#'
#' @param a1,b1 numeric: location & scale parameters for background mortality,
#' respectively
#' @param d1 character: probability distribution chosen to describe data
#' @return a vector
#' @details The expected average longevity is calculated as the integral from
#'   zero to infinity of the cumulative survival function for background
#'   mortality,
#'   given values of a1, b1, d1
#' @seealso \code{\link{av_long_infected}}
#' @examples
#'   av_long_uninfected(a1 = 3.0, b1 = 0.6, d1 = "Weibull") #17.947
av_long_uninfected <- function(a1, b1, d1 = ""){

  surv_controls <- function(t, a1, b1, d1 = ""){

    z1 <- P_get_zx(t, a1, b1, d1)
    S1 <- P_get_Sx(t, z1, d1)

    return(S1)
  }
  output <- stats::integrate(surv_controls, 0, Inf, a1, b1, d1)
  return(output)
}

#' Average longevity: estimate for infected hosts
#'
#' Calculates expected longevity of infected hosts due to background mortality
#' and mortality due to infection
#' @param a1,b1 numeric: location & scale parameters for background mortality,
#' respectively
#' @param a2,b2 numeric: location & scale parameters for mortality due to
#'  infection, respectively
#' @param d1,d2 character: probability distributions to describe background
#'  mortality and mortality due to infection, respectively
#' @return a vector
#' @details The expected average longevity is calculated as the integral
#' from zero to infinity for the product of the cumulative survival
#' functions for background mortality and mortality due to infection,
#' given values of a1, b1, d1, a2, b2, d2
#' @seealso \code{\link{av_long_uninfected}}
#' @examples
#'   av_long_infected(
#'     a1 = 3.0, b1 = 0.6, d1 = "Weibull",
#'     a2 = 2.5, b2 = 0.5, d2 = "Frechet") # 12.156
av_long_infected <- function(a1, b1, a2, b2, d1 = "", d2 = ""){

  surv_infected <- function(t, a1, b1, a2, b2, d1 = "", d2 = ""){

    z1 <- P_get_zx(t, a1, b1, d1)
    z2 <- P_get_zx(t, a2, b2, d2)

    S1 <- P_get_Sx(t, z1, d1)
    S2 <- P_get_Sx(t, z2, d2)

    S1S2 <- S1 * S2

    return(S1S2)
  }
  output <- stats::integrate(surv_infected, 0, Inf, a1, b1, a2, b2, d1, d2)
  return(output)
}






#' Expected time of death: infected hosts
#'
#' Time when infected hosts are expected to have died due to their cumulative
#' exposure to background mortality and mortality due to infection.
#'
#' It is the time t when, H1[t] + H2[t] = 1,
#' and cumulative survival is, S1[t].S2[t] = exp(-1) = 0.367
#' @param a1,b1 location & scale parameters for background mortality
#' @param a2,b2 location & scale parameters for mortality due to infection
#' @param d1,d2 character: name of probability distributions describing
#' background mortality and mortality due to infection, respectively
#' @param tmax numeric. Maximum value of t for which expected time of death
#' is searched; defaults to 100. Minimum time is zero (0)
#' @return numeric
#'@seealso \code{\link{etd_uninfected}}
#' @examples
#' print(etd_infected(a1 = 2, b1 = 0.5, a2 = 30, b2 = 5,
#'                    d1 = "Weibull", d2 = "Gumbel")) # 7.34
#' print(etd_infected(a1 = 20, b1 = 5, a2 = 3, b2 = 0.6,
#'                    d1 = "Gumbel", d2 = "Frechet")) # 17.84
#' print(etd_infected(a1 = 3, b1 = 0.6, a2 = 2, b2 = 0.5,
#'                    d1 = "Frechet", d2 = "Weibull")) # 7.37
etd_infected <- function(a1, b1, a2, b2, d1 = "", d2 = "", tmax = 100){

  t <- 1

  # get expression for H1 & replace ax, bx with a1, b1
  H1 <- P_get_func_by_name_of_distribution_and_function(d1, survfunc = "H")
  H1a <- lapply(H1, function(x) gsub("ax", "a1", x))
  H1ab <- lapply(H1a, function(x) gsub("bx", "b1", x))

  # get expression for H2 & replace ax, bx with a2, b2
  H2 <- P_get_func_by_name_of_distribution_and_function(d2, survfunc = "H")
  H2a <- lapply(H2, function(x) gsub("ax", "a2", x))
  H2ab <- lapply(H2a, function(x) gsub("bx", "b2", x))

  # add expressions H1 & H2
  H1ab_H2ab <- paste(H1ab, H2ab, sep = " + ")

  # to solve when H1 + H2 = 1, need to add text for substract '1'
  H1ab_H2ab_1 <- paste(H1ab_H2ab, "1", sep = " - ")

  # single factor expression in t to for uniroot to solve
  eqn_t <- function(t){eval(parse(text = H1ab_H2ab_1))}

  solve_t <- stats::uniroot(eqn_t, c(0, tmax))$root
}



#' Expected time of death: uninfected hosts
#'
#' Time when uninfected hosts are expected to have died due to their
#' cumulative exposure to background mortality.
#'
#' It is the time t when, H1[t] = 1, and
#' cumulative survival is, S1[t] = exp(-1) = 0.367
#' @param a1,b1 location & scale parameters for background mortality
#' @param d1 character: name of probability distribution describing
#' background mortality
#' @param tmax numeric. Maximum value of t for which expected time of death
#' is searched; defaults to 100. Minimum time is zero (0)
#' @return numeric
#'@seealso \code{\link{etd_infected}}
#' @examples
#' print(etd_uninfected(a1 = 2, b1 = 0.5, d1 = "Weibull")) # 7.38
#' print(etd_uninfected(a1 = 20, b1 = 5, d1 = "Gumbel")) # 20
#' print(etd_uninfected(a1 = 3, b1 = 0.6, d1 = "Frechet")) # 32.06
etd_uninfected <- function(a1, b1, d1 = "", tmax = 100){

  t <- 1

  # get expression for H1 & replace ax, bx with a1, b1
  H1 <- P_get_func_by_name_of_distribution_and_function(d1, survfunc = "H")
  H1a <- lapply(H1, function(x) gsub("ax", "a1", x))
  H1ab <- lapply(H1a, function(x) gsub("bx", "b1", x))

  # to solve when H1 = 1, need to add text for substract '1'
  H1ab_1 <- paste(H1ab, "1", sep = " - ")

  # single factor expression in t to for uniroot to solve
  eqn_t <- function(t){eval(parse(text = H1ab_1))}

  solve_t <- stats::uniroot(eqn_t, c(0, tmax))$root
}




#' Approximate 95\% confidence intervals for virulence
#'
#' Function calculating the 95\% confidence intervals
#' for a hazard function based on the variance and covariance
#' of its location and scale parameters.
#'
#' The approach is based on the interval being estimated as a complementary
#' log-log function of the hazard function, h(t), with the variance
#' of virulence being estimated by the Delta method applied to log(h[t]).
#'
#' @param a2 numeric. Estimated value of location parameter describing
#' mortality due to infection
#' @param b2 numeric. Estimated value of scale parameter describing
#' mortality due to infection
#' @param var_a2 numeric. Estimated variance of location parameter describing
#' mortality due to infection
#' @param var_b2 numeric. Estimated variance of scale parameter describing
#' mortality due to infection
#' @param cov_a2b2 numeric. Estimated covariance of location and scale
#' parameters above
#' @param d2 character. Probability distribution assumed to describe
#'   virulence; Weibull, Gumbel or Fréchet
#' @param tmax maximum time virulence will be calculated for. Default value;
#'   tmax = 21
#' @return matrix containing estimates of virulence over time ± approx. 95\%
#'   confidence intervals
#'
#' @examples
#'
#' # the values, variance and covariance of the location and scale parameters
#' # [a2,a2] describing mortality due to infection were estimated as;
#' # a2 = 2.5807642
#' # b2 = 0.1831328
#' # var_a2 = 0.0008196927
#' # var_b2 = 0.0010007282
#' # cov_a2b2 = -0.0003119921
#'
#'  ci_matrix01 <- conf_ints_virulence(
#'    a2 = 2.5807642,
#'    b2 = 0.1831328,
#'    var_a2 = 0.0008196927,
#'    var_b2 = 0.0010007282,
#'    cov_a2b2 = -0.0003119921,
#'    d2 = "Weibull",
#'    tmax = 15)
#'
#'  tail(ci_matrix01)
#'
#'  plot(ci_matrix01[, 't'], ci_matrix01[, 'h2'],
#'    type = 'l', col = 'red',
#'    xlab = 'time', ylab = 'virulence (± 95% ci)')
#'    lines(ci_matrix01[, 'lower_ci'], col = 'grey')
#'    lines(ci_matrix01[, 'upper_ci'], col = 'grey')
#'
conf_ints_virulence <- function(
  a2 = a2,
  b2 = b2,
  var_a2 = var_a2,
  var_b2 = var_b2,
  cov_a2b2 = cov_a2b2,
  d2 = "",
  tmax = 21){

  # a2, b2 are location and scale parameters describing
  # mortality due to infection, respectively

  t <- 1:tmax

  #  a2 <- mle2object@coef[[3]]
  #  b2 <- mle2object@coef[[4]]

  # var_a2 <- mle2object@vcov[3, 3]
  # var_b2 <- mle2object@vcov[4, 4]
  # cov_a2b2 <- mle2object@vcov[3, 4]

  a2 <- a2
  b2 <- b2
  var_a2 <- var_a2
  var_b2 <- var_b2
  cov_a2b2 <- cov_a2b2

  # get expression for hazard function
  hazard_function <- P_hazard_expression(d2)

  # calculate partial derivatives of 'h2' by 'a2' and 'b2'
  dh2_dt <- stats::deriv(hazard_function, c('a2', 'b2'))

  h2 <- eval(dh2_dt)

  dh2 <- attr(h2, 'gradient')

  colnames(dh2) <- c('dh2_da2', 'dh2_db2')

  dh2 <- cbind(t, h2, dh2)

  # calculate variance of hazard function
  var_h2 <- dh2[,'dh2_da2']^2 * var_a2 +
    dh2[,'dh2_db2']^2 * var_b2 +
    2 * dh2[,'dh2_da2'] * dh2[,'dh2_db2'] * cov_a2b2

  sd_h2 <- sqrt(var_h2)

  dh2 <- cbind(dh2, sd_h2)

  lower_ci <- dh2[,'h2'] * exp(-1.96 * dh2[,'sd_h2'] / dh2[,'h2'])
  upper_ci <- dh2[,'h2'] * exp( 1.96 * dh2[,'sd_h2'] / dh2[,'h2'])

  dh2 <- cbind(dh2, lower_ci, upper_ci)

}



#' Checks data are correctly described for models
#'
#' Function checking 'time', 'censor' and 'infected treatment'
#' columns in data.frame are correctly specified for the
#' likelihood functions in this package.
#'
#' An error message is also given if data contain rows 'NA';
#' these need removing.
#'
#' NB error messages triggered on 1st encounter with a fault.
#' Correct and check revised data for other faults.
#' @param data Name of data.frame to be checked
#' @param time Name of column with time of event data (death or censoring);
#'
#' requires time > 0 and numeric.
#' @param censor Name of column containing censor status;
#'
#'  '1' censored,
#'
#'  '0' uncensored or died during experiment
#'
#'  (needs to be numeric)
#' @param infected_treatment Name of column for infection treatment status;
#'
#' '1' a treatment exposed to infection,
#'
#' '0' treatment not exposed to infection
#'
#' (needs to be numeric).
#'
#' NB Different likelihood models make different assumptions as to whether
#' all the individuals in an infected treatment are infected
#' or not (c.f. \code{\link{nll_exposed_infected}})
#' @return Message, 'Checks completed' or error message
#' @examples
#' \donttest{
#' # view head of data.frame to be checked
#'   head(data_blanford, 3)
#'
#' # specify data.frame and names of columns for;
#'   # time, censor status, infection status
#'   check_data(data = data_blanford,
#'              time = t,
#'              censor = censor,
#'              infected_treatment = inf)
#'
#' # create data.frame 't_zero' with t = 0 in first row of data.frame
#' t_zero <- data_blanford
#' t_zero[1, 8] <- 0
#' head(t_zero, 3)
#'
#' check_data(data = t_zero,
#'            time = t,
#'            censor = censor,
#'            infected_treatment = inf)
#'
#' # correct '0' and make first row of column 'censor' = NA
#'   t_zero[1, 8] <- 1 ; t_zero[1, 5] <- NA ; head(t_zero, 3)
#'
#' check_data(data = t_zero,
#'            time = t,
#'            censor = censor,
#'            infected_treatment = inf)
#'
#' # remove row(s) with 'NA'
#'   t_zero_II <- na.omit(t_zero) # NB applies to whole data.frame
#'   check_data(data = t_zero_II,
#'            time = t,
#'            censor = censor,
#'            infected_treatment = inf)
#'
#' }
check_data <- function(data = data,
                       time = time,
                       censor = censor,
                       infected_treatment = infected_treatment){

  time <- data[[deparse(substitute(time))]]
  censor <- data[[deparse(substitute(censor))]]
  infected_treatment <- data[[deparse(substitute(infected_treatment))]]

  stopifnot(is.numeric(time),
            is.numeric(censor),
            is.numeric(infected_treatment),
            !is.na(time),
            !is.na(censor),
            !is.na(infected_treatment),
            ((censor == 0) | (censor == 1)), # censor = (0 or 1)
            ((infected_treatment == 0) | (infected_treatment == 1)),
            time > 0 # time > 0 ; avoids problems with log(t)
  )
  print('Checks completed')
}


