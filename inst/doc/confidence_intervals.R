## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(anovir)

## ----include = FALSE----------------------------------------------------------

  sy_times <- c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,23,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,7,21,22,23,7,14,21,22,23)
  sy_censor <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
  sy_inf <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1)
  sy_fq <- c(5,8,8,4,3,2,4,3,3,1,3,3,4,2,3,4,5,2,2,1,7,4,8,9,4,5,4,2,3,6,10,7,13,5,11,5,9,11,5,5,3,4,2,3,4,205,143,6,7,4,103,66)

  data_sy <- data.frame(sy_times, sy_inf, sy_censor, sy_fq)
  
  data_sy$fq <- data_sy$sy_fq

## ----include = FALSE, eval = TRUE, echo = FALSE, message = FALSE, warning = FALSE----

    m01_prep_function <- function(a1, b1, a2, b2){
      nll_basic(
        a1, b1, a2, b2,
        data = data_sy,
        time = sy_times,
        censor = sy_censor,
        infected_treatment = sy_inf,
        d1 = 'Weibull', d2 = 'Weibull')
        }


## ----include = FALSE, eval = TRUE, echo = FALSE, message = FALSE, warning = FALSE----
    m02 <- mle2(m01_prep_function,
             start = list(a1 = 4.8, a2 = 3.6, b2 = 0.6),
             fixed = list(b1 = 1.0)
             )

## ----include = TRUE-----------------------------------------------------------
  coef(m02)

## ----include = TRUE-----------------------------------------------------------
  vcov(m02)

## ----include = TRUE-----------------------------------------------------------
  mu_a1 <- coef(m02)[[1]] ; mu_a1
  var_a1 <- vcov(m02)[1,1] ; var_a1

  h1t <- 1/(exp(mu_a1)) ; h1t
  var_h1t <- (-1/exp(mu_a1))^2 * var_a1 ; var_h1t

  lower_ci <- h1t - 1.96*sqrt(var_h1t)
  upper_ci <- h1t + 1.96*sqrt(var_h1t)
  lower_ci ; h1t ; upper_ci

## ----include = TRUE-----------------------------------------------------------
  lower_ci2 <- h1t*exp(-1.96*sqrt(var_h1t)/h1t)
  upper_ci2 <- h1t*exp(+1.96*sqrt(var_h1t)/h1t)

  lower_ci2 ; h1t ; upper_ci2

## ----include = TRUE, warning = FALSE------------------------------------------

library(anovir)

data01 <- data_blanford 

data01 <- subset(data01,
    (data01$block == 5) & (
      (data01$treatment == 'cont') |
      (data01$treatment == 'Ma07') 
    ) &
    (data01$day > 0)
    )

# create column 'g' as index of infected treatment
  data01$g <- data01$inf

  m01_prep_function <- function(a1, b1, a2, b2){
    nll_basic(a1, b1, a2, b2, data = data01, time = t, censor = censor,
              infected_treatment = g, d1 = 'Weibull', d2 = 'Fréchet')
    }

  m01 <- mle2(m01_prep_function,
              start = list(a1 = 2, b1 = 0.5, a2 = 3, b2 = 1)
              )

  summary(m01)


## ----includes = TRUE----------------------------------------------------------

# view estimates
  coef(m01) ; vcov(m01)

# assign means, variances & covariances 
# means
  a2 <- m01@coef[[3]]
  b2 <- m01@coef[[4]]
# variances
  var_a2 <- m01@vcov[3,3]
  var_b2 <- m01@vcov[4,4]
# covariances
  cov_a2_b2 <- m01@vcov[3,4]

# specify timespan to calculate
  t <- 1:14

# write an 'expression' for the Fréchet hazard function
  # in terms of 'a2', 'b2' 
  hazard_expression <- expression(
    (1 / (b2 * t)) * exp(-((log(t) - a2) / b2) - exp(-((log(t) - a2) / b2))) / (1 - exp(-exp(-((log(t) - a2) / b2))))
  )

# prepare expression for 'stats::deriv` 
  # NB needs names of variables for partial differentiation 
  dh2_dt <- stats::deriv(hazard_expression, c('a2', 'b2'))

  dh2_dt
  
# evaluate (calculate) hazard function 'h2'
  h2 <- eval(dh2_dt)

# collect estimate of hazard function 'h2' 
# and 'gradients' of partial derivatives
  dh2 <- attr(h2, 'gradient')

# rename columns for clarity
  colnames(dh2) <- c('dh2_da2', 'dh2_db2')

# collect time, estimates of hazard function at time 't'
# and partial derivatives 
  dh2 <- cbind(t, h2, dh2)  

# calculate variance of hazard function
# using delta function  
  var_h2 <- 
    dh2[,'dh2_da2']^2 * var_a2 +
    dh2[,'dh2_db2']^2 * var_b2 +
    2 * dh2[,'dh2_da2'] * dh2[,'dh2_db2'] * cov_a2_b2

# calculate standard deviation of estimate
  sd_h2 <- sqrt(var_h2)

# bind to other estimates
  dh2 <- cbind(dh2, sd_h2)

# calculate lower & upper bounds of 95% confidence interval
  lower_ci <- dh2[,'h2'] * exp(-1.96 * dh2[,'sd_h2'] / dh2[,'h2'])
  upper_ci <- dh2[,'h2'] * exp( 1.96 * dh2[,'sd_h2'] / dh2[,'h2'])

# bind together
  dh2_ma07 <- cbind(dh2, lower_ci, upper_ci)


## ----includes = TRUE, echo = FALSE, fig.height = 3.25, fig.width = 3.25, fig.show = 'hold'----

# plot data
  plot(dh2_ma07[,'t'], dh2_ma07[,'h2'],
    type = 'l', col = 'red', lty = 'solid',
    xlim = c(0, 14), ylim = c(0, 0.6),
    xlab = 'Time(days)', ylab = 'h(t) ±95% c.i.',
    main = 'Virulence: Ma07',
    axes = FALSE
    )
  axis(side = 1, seq(0, 14, by = 7))
  axis(side = 2, seq(0, 0.6, by = 0.1))
  lines(dh2_ma07[,'t'], dh2_ma07[,'lower_ci'],
    type = 'l', col = 'grey', lty = 'solid'
    )
  lines(dh2_ma07[,'t'], dh2_ma07[,'upper_ci'],
    type = 'l', col = 'grey', lty = 'solid'
    )

## ----include = TRUE-----------------------------------------------------------

# values for a2, b2, var_a2, var_b2, cov_a2b2 were assigned 
# from results of 'bbmle::mle2' above
ls()

# tail of calculations above
tail(dh2_ma07)

# specify mle2object 'm01', Frechet distribution & maximum time
ci_matrix01 <- conf_ints_virulence(
  a2 = a2,
  b2 = b2,
  var_a2 = var_a2,
  var_b2 = var_b2,
  cov_a2b2 = cov_a2_b2,
  d2 = "Frechet", 
  tmax = 14)  

# tail of calculations from 'conf_ints_virulence'
tail(ci_matrix01)

# are the results the same? 
identical(dh2_ma07, ci_matrix01)


## ----includes = TRUE, warning = FALSE-----------------------------------------

library(anovir)

# get & rename data
  data01 <- data_lorenz 
  head(data01)

# create new version of function 'nll_basic'
  nll_basic2 <- nll_basic

# check location/definition of parameter function a2 ('pfa2')
  body(nll_basic2)[[4]]

# replace 'a2' making 'pfa2' a linear function of log(dose)
  body(nll_basic2)[[4]] <- substitute(pfa2 <- a2i + a2ii * log(data01$Infectious.dose))

# check new version of 'pfa2'
  body(nll_basic2)[[4]]

# update formals, including names for time, censor, etc..
  formals(nll_basic2) <- alist(a1 = a1, b1 = b1, a2i = a2i, a2ii = a2ii, b2 = b2,
    data = data01, time = t, censor = censored, infected_treatment = g,
    d1 = 'Gumbel', d2 = 'Weibull')

# prepare 'prep_function' for analysis
  m01_prep_function <- function(a1, b1, a2i, a2ii, b2){
    nll_basic2(a1, b1, a2i, a2ii, b2)
  }

# send to 'bbmle::mle2' with starting values for variables
  m01 <- mle2(m01_prep_function,
      start = list(a1 = 23, b1 = 4, a2i = 4, a2ii = -0.1, b2 = 0.2)
      )

# summary of results
  summary(m01)


## ----includes = TRUE----------------------------------------------------------

# collect & assign results in mleobject 'm01'
  coef(m01)
  vcov(m01)

# means
           a2i <- m01@coef[[3]]
          a2ii <- m01@coef[[4]]
            b2 <- m01@coef[[5]]
# variances
       var_a2i <- m01@vcov[3,3]
      var_a2ii <- m01@vcov[4,4]
        var_b2 <- m01@vcov[5,5]
# covariances
  cov_a2i_a2ii <- m01@vcov[3,4]
    cov_a2i_b2 <- m01@vcov[3,5]
   cov_a2ii_b2 <- m01@vcov[4,5]

# specify timespan to calculate
  t <- 1:28

# specify dose treatment to calculate
  dose <- 5000

# write an 'expression' for hazard function
  hazard_expression <- expression(
    1 / (b2 * t) * exp(((log(t) - (a2i + a2ii * log(dose)))) / b2)
    )

# prepare expression for 'stats::deriv` 
# has names of variables for partial differentiation 
  dh2_dt <- stats::deriv(hazard_expression, c('a2i', 'a2ii', 'b2'))

# evaluate (calculate) hazard function 'h2'
  h2 <- eval(dh2_dt)

# collect estimate of hazard function 'h2' 
# and 'gradients' of partial derivatives
  dh2 <- attr(h2, 'gradient')

# rename columns for clarity
  colnames(dh2) <- c('dh2_da2i', 'dh2_da2ii', 'dh2_db2')

# collect time, estimates of hazard function at time 't'
# and partial derivatives 
  dh2 <- cbind(t, h2, dh2)

# calculate variance of hazard function
# using delta function  
  var_h2 <- 
    dh2[,'dh2_da2i']^2 * var_a2i +
    dh2[,'dh2_da2ii']^2 * var_a2ii +
    dh2[,'dh2_db2']^2 * var_b2 +
    2 * dh2[,'dh2_da2i'] * dh2[,'dh2_da2ii'] * cov_a2i_a2ii +
    2 * dh2[,'dh2_da2i'] * dh2[,'dh2_db2'] * cov_a2i_b2 +
    2 * dh2[,'dh2_da2ii'] * dh2[,'dh2_db2'] * cov_a2ii_b2 

# calculate standard deviation of estimate
  sd_h2 <- sqrt(var_h2)

# bind to other estimates
  dh2 <- cbind(dh2, sd_h2)

# calculate lower & upper bounds of 95% confidence interval
  lower_ci <- dh2[,'h2'] * exp(-1.96 * dh2[,'sd_h2'] / dh2[,'h2'])
  upper_ci <- dh2[,'h2'] * exp( 1.96 * dh2[,'sd_h2'] / dh2[,'h2'])

# bind together
  dh2_5000 <- cbind(dh2, lower_ci, upper_ci)


## ----includes = TRUE, echo = FALSE--------------------------------------------

# specify dose treatment to calculate
  dose <- 160000

# write an 'expression' for hazard function
  hazard_expression <- expression(
    1 / (b2 * t) * exp(((log(t) - (a2i + a2ii * log(dose)))) / b2)
    )

# prepare expression for 'stats::deriv` 
# has names of variables for partial differentiation 
  dh2_dt <- stats::deriv(hazard_expression, c('a2i', 'a2ii', 'b2'))

# evaluate (calculate) hazard function 'h2'
  h2 <- eval(dh2_dt)

# collect estimate of hazard function 'h2' 
# and 'gradients' of partial derivatives
  dh2 <- attr(h2, 'gradient')

# rename columns for clarity
  colnames(dh2) <- c('dh2_da2i', 'dh2_da2ii', 'dh2_db2')

# collect time, estimates of hazard function at time 't'
# and partial derivatives 
  dh2 <- cbind(t, h2, dh2)

# calculate variance of hazard function
# using delta function  
  var_h2 <- 
    var_a2i * dh2[,'dh2_da2i']^2 +
    var_a2ii * dh2[,'dh2_da2ii']^2 +
    var_b2 * dh2[,'dh2_db2']^2 +
    2 * cov_a2i_a2ii * dh2[,'dh2_da2i'] * dh2[,'dh2_da2ii'] +
    2 * cov_a2i_b2 * dh2[,'dh2_da2i'] * dh2[,'dh2_db2'] +
    2 * cov_a2ii_b2 * dh2[,'dh2_da2ii'] * dh2[,'dh2_db2']

# calculate standard deviation of estimate
  sd_h2 <- sqrt(var_h2)

# bind to other estimates
  dh2 <- cbind(dh2, sd_h2)

# calculate lower & upper bounds of 95% confidence interval
  lower_ci <- dh2[,'h2'] * exp(-1.96 * dh2[,'sd_h2'] / dh2[,'h2'])
  upper_ci <- dh2[,'h2'] * exp( 1.96 * dh2[,'sd_h2'] / dh2[,'h2'])

# bind together
  dh2_160000 <- cbind(dh2, lower_ci, upper_ci)

## ----includes = TRUE, echo = FALSE, fig.height = 3.25, fig.width = 3.25, fig.show = 'hold'----

# plot data
  plot(dh2_5000[,'t'], dh2_5000[,'h2'],
    type = 'l', col = 'red', lty = 'solid',
    xlim = c(0, 28), ylim = c(0, 1.2),
    xlab = 'Time(days)', ylab = 'h(t) ±95% c.i.',
    main = 'Virulence: dose 5000',
    axes = FALSE
    )
  axis(side = 1, seq(0, 28, by = 7))
  axis(side = 2, seq(0, 1.2, by = 0.2))
  lines(dh2_5000[,'t'], dh2_5000[,'lower_ci'],
    type = 'l', col = 'grey', lty = 'solid'
    )
  lines(dh2_5000[,'t'], dh2_5000[,'upper_ci'],
    type = 'l', col = 'grey', lty = 'solid'
    )

# plot data
  plot(dh2_160000[,'t'], dh2_160000[,'h2'],
    type = 'l', col = 'red', lty = 'solid',
    xlim = c(0, 28), ylim = c(0, 1.2),
    xlab = 'Time(days)', ylab = 'h(t) ±95% c.i.',
    main = 'Virulence: dose 160000',
    axes = FALSE
    )
  axis(side = 1, seq(0, 28, by = 7))
  axis(side = 2, seq(0, 1.2, by = 0.2))
  lines(dh2_160000[,'t'], dh2_160000[,'lower_ci'],
    type = 'l', col = 'grey', lty = 'solid'
    )
  lines(dh2_160000[,'t'], dh2_160000[,'upper_ci'],
    type = 'l', col = 'grey', lty = 'solid'
    )


## ----includes = TRUE----------------------------------------------------------

library(anovir)

  data01 <- subset(data_blanford,
    (data_blanford$block == 3) &
      ((data_blanford$treatment == 'cont') |
      (data_blanford$treatment == 'Bb06')) &
    (data_blanford$day > 0)
    )

    l01_prep_function_log <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2){
      nll_basic_logscale(
        a1 = a1, b1 = b1, a2 = a2, b2 = b2,
        data = data01,
        time = t,
        censor = censor,
        infected_treatment = inf,
        d1 = 'Weibull', d2 = 'Weibull')
        }

    l01 <- mle2(l01_prep_function_log,
             start = list(a1 = 1, b1 = 1, a2 = 1, b2 = 1)
             )

    summary(l01)

## ----include = TRUE-----------------------------------------------------------

# collect and assign results
  coef(l01)
  vcov(l01)

        a2 <- coef(l01)[[3]]
        b2 <- coef(l01)[[4]]

    var_a2 <- vcov(l01)[3,3]
    var_b2 <- vcov(l01)[4,4]

  cov_a2_b2 <- vcov(l01)[3,4]

# set timescale for calculations
  t <- 1:15

# NB need to exponentiate the terms with 'a2' and 'b2'
# in hazard expression
  hazard_expression <- expression(
    (1/(exp(b2) * t)) * exp((log(t) - exp(a2)) / exp(b2))
  )

# remaining steps are the same as in Example 2 above
  dh2_dt <- stats::deriv(hazard_expression, c('a2', 'b2'))
  h2 <- eval(dh2_dt)
  dh2 <- attr(h2, 'gradient')
  colnames(dh2) <- c('dh2_da2', 'dh2_db2')
  dh2 <- cbind(t, h2, dh2) 
  var_h2 <- 
    dh2[,'dh2_da2']^2 * var_a2 +
    dh2[,'dh2_db2']^2 * var_b2 +
    2 * dh2[,'dh2_da2'] * dh2[,'dh2_db2'] * cov_a2_b2
  sd_h2 <- sqrt(var_h2)
  dh2 <- cbind(dh2, sd_h2)
  lower_ci <- dh2[,'h2'] * exp(-1.96 * dh2[,'sd_h2'] / dh2[,'h2'])
  upper_ci <- dh2[,'h2'] * exp( 1.96 * dh2[,'sd_h2'] / dh2[,'h2'])

  dh2_backlog <- cbind(dh2, lower_ci, upper_ci)
  dh2_backlog <- round(dh2_backlog,4)
  dh2_backlog

  # the confidence intervals are the same as those
  # estimated on the help page of
  # `conf_ints_virulence` for variables that are 
  # not assumed to be on a logscale


