## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(anovir)

## ----include = TRUE, warning = FALSE------------------------------------------

data01 <- data_lorenz

# create column 'ref_vir' in dataframe coded 0, 1 and 2
# for control, reference dose, and other dose treatments, respectively

data01$ref_vir <- ifelse(data01$g == 0, 0,
  ifelse(data01$g == 1, ifelse(data01$Infectious.dose == 5000, 1, 2), 0)
  )

# copy and modify 'nll_proportional_virulence' function
# where the virulence in the reference treatment is 'theta'
# which is multiplied by 'th10', 'th20', etc... in the other dose treatements

nll_proportional_virulence2 <- nll_proportional_virulence
  body(nll_proportional_virulence2)[[6]] <- substitute(pfth <- theta *
  ifelse(data01$Infectious.dose == 10000, th10,
  ifelse(data01$Infectious.dose == 20000, th20,
  ifelse(data01$Infectious.dose == 40000, th40,
  ifelse(data01$Infectious.dose == 80000, th80,
  ifelse(data01$Infectious.dose == 160000, th160,
  exp(0)
  )))))
  )

# update formals 

formals(nll_proportional_virulence2) <- alist(
  a1 = a1,  b1 = b1,  a2 = a2,  b2 = b2,
  theta = theta,
  th10 = th10, th20 = th20, th40 = th40, th80 = th80, th160 = th160,
  data = data01,
  time = t,
  censor = censored,
  infected_treatment = g,
  reference_virulence = ref_vir,
  d1 = 'Gumbel', d2 = 'Weibull')

# write 'prep function'

m01_prep_function <- function(
  a1, b1, a2, b2, theta, th10, th20, th40, th80, th160){
  nll_proportional_virulence2(
  a1, b1, a2, b2, theta, th10, th20, th40, th80, th160)
  }

# send 'prep function' to mle2
# NB fixing 'theta = 1' scales virulence in other treatments relative to
# that in the reference treatment

m01 <- mle2(m01_prep_function,
  start = list(
    a1 = 24, b1 = 5, a2 = 4, b2 = 0.2,
    theta = 1,
    th10 = 1, th20 = 1, th40 = 1, th80 = 1, th160 = 1),
    fixed = list(theta = 1)
    )

summary(m01)


## ----include = TRUE-----------------------------------------------------------
data01 <- data_lorenz

head(data01)


nll_basic_logscale2 <- nll_basic_logscale 

body(nll_basic_logscale2)[[4]]


## ----include = TRUE-----------------------------------------------------------

body(nll_basic_logscale2)[[4]] <- substitute(
  pfa2 <- exp(a2i) - exp(a2ii)*log(data01$Infectious.dose)
  )

body(nll_basic_logscale2)[[4]]

formals(nll_basic_logscale2) <- alist(
  a1 = a1, b1 = b1, a2i = a2i, a2ii = a2ii, b2 = b2,
  data = data01,
  time = t,
  censor = censored,
  infected_treatment = g,
  d1 = 'Gumbel', d2 = 'Weibull'
)

m01_prep_function <- function(a1, b1, a2i, a2ii, b2){
  nll_basic_logscale2(a1, b1, a2i, a2ii, b2)
  }

m01 <- mle2(m01_prep_function,
  start = list(a1 = 3, b1 = 1.5, a2i = 1.4, a2ii = -3, b2 = -1.5)
  )

summary(m01)

exp(coef(m01))


## ----include = TRUE-----------------------------------------------------------

nll_basic2 <- nll_basic 

body(nll_basic2)[[4]]

body(nll_basic2)[[4]] <- substitute(
  pfa2 <- a2i + a2ii*log(data01$Infectious.dose)
  )

body(nll_basic2)[[4]]

formals(nll_basic2) <- alist(
  a1 = a1, b1 = b1, a2i = a2i, a2ii = a2ii, b2 = b2,
  data = data01,
  time = t,
  censor = censored,
  infected_treatment = g,
  d1 = 'Gumbel', d2 = 'Weibull'
)


m02_prep_function <- function(a1, b1, a2i, a2ii, b2){
  nll_basic2(a1, b1, a2i, a2ii, b2)
  }

m02 <- mle2(m02_prep_function,
  start = list(a1 = 23, b1 = 4.5, a2i = 3.8, a2ii = -0.1, b2 = 0.4)
  )

summary(m02)



