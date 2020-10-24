
#' Negative log-likelihood function: recovery model, no background mortality
#'
#' Function returning negative log-likelihood (nll) for patterns of survival in
#' infected and uninfected treatments, when infected hosts can recover from
#' infection and there is no background mortality.
#'
#' This model assumes all the hosts in an infected treatment are all initially
#' infected, and they can all potentially recover from infection. Uninfected,
#' infected and recovered hosts are assumed to not experience any background
#' mortality during the period of the experiment. No assumptions are made as to
#' whether recovered hosts are still infected or infectious. It is also assumed
#' that the timing of recovery from infection is not directly observed, but and
#' individual's infected/recovery status can be determined after they have died
#' or been censored.
#'
#' The probability that an infection 'survives' over time, i.e., the host does
#' not recover from infection, is assumed to follow a probability distribution
#' which acts independently of the probability distributions determining
#' background mortality or mortality due to infection.
#'
#' This function only estimates location and scale parameters as constants, it
#' is not designed to estimate them as functions.
#'
#' The nll function also requires the data to be specified in a specific format
#' and requires columns detailing when and how many control individuals were
#' right-censored, even though these individuals do not contribute to the nll
#' estimated; see vignettes for details.
#'
#'
#' @param a2,b2 location & scale parameters for mortality due to infection
#' @param a3,b3 location & scale parameters for how long infection 'survives'
#' @param data a data.frame with the data
#' @param d2,d3 probability distributions for background mortality, mortality
#'   due to infection & how long infection 'survives' ("Weibull", "Gumbel",
#'   "Frechet")
#' @return numeric
#' @section Warning: requires the data to be specified in a specific format;
#' see vignette 'data format' for details
#' @examples
#' \donttest{
#' # NB the data to analyse needs to be in a data frame of a specific form
#'      head(recovery_data_II)
#'
#' # step #1: prepare nll function for analysis
#'     m01_prep_function <- function(a2, b2, a3, b3){
#'       nll_recovery_II(a2, b2, a3, b3,
#'                    data = recovery_data_II, # data_recovery_II,
#'                    d2 = "Weibull", d3 = "Weibull"
#'                   )}
#'
#' # step #2: send 'prep_function' to mle2 for maximum likelihood estimation,
#'   # specifying starting values
#'     m01 <- mle2(m01_prep_function,
#'              start = list(a2 = 2, b2 = 0.5, a3 = 2, b3 = 0.5)
#'              )
#'
#'     summary(m01)
#'
#' # values used to simulate data were for the Weibull distribution;
#'   # a2 = 2.2, b2 = 0.35, a3 = 2.35, b3 = 0.35
#' }
nll_recovery_II <- function(a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                            data = data, d2 = "", d3 = ""){
                              nll <- P_recovery_II_calc_likelihood_matrix(
                              a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                              data = data, d2 = "", d3 = "")
                              }


###

P_recovery_II_calc_survival_functions <- function(
  a2 = a2, b2 = b2, a3 = a3, b3 = b3, data = data, d2 = "", d3 = ""
  ){

  # function returns matrix with values for various survival functions
  # from t[1] -> t[max] given values of a1, b1, etc.

  tmax <- max(data$t)

  calc.matrix <- matrix(0,tmax,13)

  colnames(calc.matrix) <- c("t", "f1", "f2", "f3", "S1", "S2", "S3", "h1",
                             "h2", "h3", "f1S2S3", "f2S1S3", "f3S1S2")

  for (t in 1:tmax){

    #    z1 <- P_get_zx(t, a1, b1, d1)
    z2 <- P_get_zx(t, a2, b2, d2)
    z3 <- P_get_zx(t, a3, b3, d3)

    f1 <- 0
    f2 <- P_get_fx(t, z2, b2, d2)
    f3 <- P_get_fx(t, z3, b3, d3)

    S1 <- 1
    S2 <- P_get_Sx(t, z2, d2)
    S3 <- P_get_Sx(t, z3, d3)

    h1 <- 0
    h2 <- P_get_hx(t, z2, b2, d2)
    h3 <- P_get_hx(t, z3, b3, d3)

    f1S2S3 <- 0
    f2S1S3 <- f2 * S1 * S3
    f3S1S2 <- f3 * S1 * S2

    calc.matrix[t,1] <- t
    calc.matrix[t,2] <- f1
    calc.matrix[t,3] <- f2
    calc.matrix[t,4] <- f3
    calc.matrix[t,5] <- S1
    calc.matrix[t,6] <- S2
    calc.matrix[t,7] <- S3
    calc.matrix[t,8] <- h1
    calc.matrix[t,9] <- h2
    calc.matrix[t,10] <- h3
    calc.matrix[t,11] <- f1S2S3
    calc.matrix[t,12] <- f2S1S3
    calc.matrix[t,13] <- f3S1S2
  }
  calc.matrix
}

###

P_recovery_II_calc_f3S1S2 <- function(a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                                    data = data, d2 = "", d3 = ""){

  # returns matrix with calculations for probabilities of recovery
  # while host still alive at different times given values of a1, b1, etc...

  tmax <- max(data$t)

  matrix01 <- matrix(0,tmax,tmax)
  matrix02 <- matrix(0,tmax,tmax)

  for (t in 1:tmax){
    for (u in 1:t){

      matrix01[t,u] <- u

      #      zu1 <- P_get_zu(u, a1, b1, d1)
      zu2 <- P_get_zu(u, a2, b2, d2)
      zu3 <- P_get_zu(u, a3, b3, d3)

      Su1 <- 1
      Su2 <- P_get_Su(u, zu2, d2)

      fu3 <- P_get_fu(u, zu3, b3, d3)

      matrix02[t,u] <- fu3 * Su1 * Su2
    }
  }
  matrix02
}

###

P_recovery_II_calc_St_Sr <- function(a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                                   data = data, d2 = "", d3 = ""){

  # returns matrix with calculations for probability of surviving from
  # time of recovery, S[r], to time t, S[t]: S[t]/S[r]
  # given values of a1, b1, etc.

  tmax <- max(data$t)

  matrix03 <- matrix(0,tmax,tmax)
  matrix04 <- matrix(0,tmax,tmax)

  for (t in 1:tmax){
    for (u in 1:t){
      matrix03[t,u] <- t

      #      z1t <- P_get_zx(t, a1, b1, d1)
      #      z1u <- P_get_zu(u, a1, b1, d1)

      S1t <- 1
      S1u <- 1

      matrix04[t,u] <- S1t / S1u
    }
  }
  matrix04
}

###

P_recovery_II_calc_f3S1S2_St_Sr <- function(a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                                          data = data, d2 = "", d3 = ""){

  # calculates probabilities of recovering at different times
  # and surviving to time t
  # from product of recovery_II_calc_f3S1S2 * recovery_II_St.Sr
  # given values of a1, b1, etc.

  f3S1S2 <- P_recovery_II_calc_f3S1S2(a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                                    data = data, d2 = "", d3 = "")
  St.Sr <- P_recovery_II_calc_St_Sr(a2 = a2, b2 = b2, a3 = a3, b3 = b3,
                                  data = data, d2 = "", d3 = "")

  f3S1S2_St_Sr <- f3S1S2 * St.Sr
}

###

P_recovery_II_calc_sum_f3S1S2_St_Sr <- function(a2 = a2, b2 = b2,
                                              a3 = a3, b3 = b3,
                                              data = data, d2 = "", d3 = ""){

  # calculates sum of probabilities of recovering at different times
  # and surviving to time t
  # given values of a1, b1, etc.


  tmax <- max(data$t)

  f3S1S2.St.Sr <- P_recovery_II_calc_f3S1S2_St_Sr(a2 = a2, b2 = b2,
                                                a3 = a3, b3 = b3,
                                                data = data, d2 = "", d3 = "")

  vector.sum.f3S1S2.St.Sr <- rowSums(f3S1S2.St.Sr)

  sum_f3S1S2_St_Sr <- matrix(vector.sum.f3S1S2.St.Sr, , 1)

}

###


P_recovery_II_calc_likelihood_matrix <- function(a2 = a2, b2 = b2,
                                               a3 = a3, b3 = b3,
                                               data = data, d2 = "", d3 = ""){

  # calculates likelihoods using 'recovery' sub-functions
  # values a1, b1, etc...

  # define maximum time
  tmax <- max(data$t)

  # convert data.frame(data) to matrix
  matrix01 <- as.matrix(data)

  # calculate various survival functions (defined elsewhere) based on values a1,a2,...

  survival.functions <- P_recovery_II_calc_survival_functions(
    a2 = a2, b2 = b2, a3 = a3, b3 = b3, data = data, d2 = "", d3 = "")
  f3S1S2 <- P_recovery_II_calc_f3S1S2(
    a2 = a2, b2 = b2, a3 = a3, b3 = b3, data = data, d2 = "", d3 = "")
  St.Sr <- P_recovery_II_calc_St_Sr(
    a2 = a2, b2 = b2, a3 = a3, b3 = b3, data = data, d2 = "", d3 = "")
  f3S1S2.St.Sr <- P_recovery_II_calc_f3S1S2_St_Sr(
    a2 = a2, b2 = b2, a3 = a3, b3 = b3, data = data, d2 = "", d3 = "")
  sum.f3S1S2.St.Sr <- P_recovery_II_calc_sum_f3S1S2_St_Sr(
    a2 = a2, b2 = b2, a3 = a3, b3 = b3, data = data, d2 = "", d3 = "")
  h1 <- 0

  # create & fill matrix for likelihood calculations
  # NB matrix column 13 has been modified to return '0' if column 12 == 0 ;
  # column13 = log(column 12)

  likelihood.matrix <-  matrix(0,nrow(matrix01), 15)

  likelihood.matrix[, 01] <- matrix01[, "t"]
  likelihood.matrix[, 02] <- matrix01[, 'control.d'] * 0
  likelihood.matrix[, 03] <- matrix01[, 'control.c'] * 1
  likelihood.matrix[, 04] <- matrix01[, 'infected.d'] *
    (survival.functions[, 'f1S2S3'] + survival.functions[, 'f2S1S3'])
  likelihood.matrix[, 05] <- matrix01[, 'infected.c'] *
    (survival.functions[, 'S1'] *
       survival.functions[, 'S2'] *
       survival.functions[, 'S3'])
  likelihood.matrix[, 06] <- sum.f3S1S2.St.Sr * h1 # recovered, died
  likelihood.matrix[, 07] <- matrix01[, 'recovered.d']
  likelihood.matrix[, 08] <- likelihood.matrix[, 6] * likelihood.matrix[, 7]
  likelihood.matrix[, 09] <- sum.f3S1S2.St.Sr # recovered, censored
  likelihood.matrix[, 10] <- matrix01[, 'recovered.c']
  likelihood.matrix[, 11] <- likelihood.matrix[, 9] * likelihood.matrix[, 10]
  likelihood.matrix[, 12] <- likelihood.matrix[, 2] +
    likelihood.matrix[, 3] +
    likelihood.matrix[, 4] +
    likelihood.matrix[, 5] +
    likelihood.matrix[, 11]
  likelihood.matrix[, 13] <- ifelse(likelihood.matrix[, 12] == 0, 0,
                                    log(likelihood.matrix[, 12]))
  likelihood.matrix[, 14] <- likelihood.matrix[, 13] *  matrix01[, 'fq']
  likelihood.matrix[, 15] <- sum(likelihood.matrix[, 14])

  # return negative log-likelihood calculated in likelihood matrix
  neg.log.likelihood <- -(likelihood.matrix[1, 15])

  #  separately uncheck the lines below & re-run function to see the
  #  calculations for each component contributing to likelihood calculation
  #  likelihood.matrix
  #  matrix01
  #  survival.functions
}

