## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(anovir)

## ----warning = FALSE----------------------------------------------------------
data01 <- subset(data_blanford,
  (data_blanford$block == 3) & 
  ((data_blanford$treatment == 'cont') | (data_blanford$treatment == 'Bb06')) &
  (data_blanford$day > 0)
  )

head(data01, 3)

m01_prep_function <- function(a1, b1, a2, b2){
  nll_basic(a1, b1, a2, b2,
    data = data01,
    time = t, 
    censor = censor,
    infected_treatment = inf,
    d1 = 'Weibull', d2 = 'Weibull')
  }

# starting values taken from linear regression of complementary
# log-log transformed cumulative survival data

m01 <- mle2(m01_prep_function,
  start = list(a1 = 3.343, b1 = 0.792, a2 = 2.508, b2 = 0.493)
  )

summary(m01)
confint(m01)

# log-likelihood based on estimates of linear regression 

m02<- mle2(m01_prep_function,
  start = list(a1 = 3.343, b1 = 0.792, a2 = 2.508, b2 = 0.493),
  eval.only = TRUE
  )

summary(m02)

AICc(m01, m02, nobs = sum(data01$fq))


## ----warning = FALSE----------------------------------------------------------

data01 <- data_lorenz
head(data01, 3)


### Model 1

m01_prep_function <- function(a1, b1, a2, b2){
  nll_basic(a1, b1, a2, b2,
    data = data01,
    time = t,
    censor = censored, 
    infected_treatment = g,
    d1 = 'Gumbel', d2 = 'Weibull')
  }

m01 <- mle2(m01_prep_function, 
  start = list(a1 = 20, b1 = 5, a2 = 3, b2 = 1)
  )

summary(m01)

### Model 2
# copy 'nll_basic' & make each parameter function a function of food treatment

nll_basic2 <- nll_basic 

body(nll_basic2)[[2]] <- substitute(
  pfa1 <- a1 + ifelse(data01$Food == 50, a1i, -a1i)
  )
body(nll_basic2)[[3]] <- substitute(
  pfb1 <- b1 + ifelse(data01$Food == 50, b1i, -b1i)
  )
body(nll_basic2)[[4]] <- substitute(
  pfa2 <- a2 + ifelse(data01$Food == 50, a2i, -a2i)
)
body(nll_basic2)[[5]] <- substitute(
  pfb2 <- b2 + ifelse(data01$Food == 50, b2i, -b2i)
)

# update formals 

formals(nll_basic2) <- alist(
  a1 = a1, a1i = a1i, 
  b1 = b1, b1i = b1i, 
  a2 = a2, a2i = a2i, 
  b2 = b2, b2i = b2i,
  data = data01, 
  time = t,
  censor = censored,
  infected_treatment = g,
  d1 = 'Gumbel', d2 = 'Weibull')

m02_prep_function <- function(a1, a1i, b1, b1i, a2, a2i, b2, b2i){
  nll_basic2(a1, a1i, b1, b1i, a2, a2i, b2, b2i)
  }

m02 <- mle2(m02_prep_function,
  start = list(
    a1 = 23.2, a1i = 0,
    b1 = 4.6, b1i = 0,
    a2 = 3.0, a2i = 0, 
    b2 = 0.2, b2i = 0)
  )

summary(m02)

### Model 3
# create columns in data frame
# with dummary variables for dose treatments

data01$d5 <- ifelse(data01$Infectious.dose == 5000, 1, 0)
data01$d10 <- ifelse(data01$Infectious.dose == 10000, 1, 0)
data01$d20 <- ifelse(data01$Infectious.dose == 20000, 1, 0)
data01$d40 <- ifelse(data01$Infectious.dose == 40000, 1, 0)
data01$d80 <- ifelse(data01$Infectious.dose == 80000, 1, 0)
data01$d160 <- ifelse(data01$Infectious.dose == 160000, 1, 0)

head(data01)

# make parameter functions 'pfa2' and 'pfb2' functions of dose
# using columns in data01

nll_basic3 <- nll_basic 

body(nll_basic3)[[4]] <- substitute(
pfa2 <- a2 + a5 * data01$d5 
           + a10 * data01$d10 
           + a20 * data01$d20 
           + a40 * data01$d40 
           + a80 * data01$d80 
           - (a5 + a10 + a20 + a40 + a80) * data01$d160
           )

body(nll_basic3)[[5]] <- substitute(
pfb2 <- b2 + b5 * data01$d5 
           + b10 * data01$d10 
           + b20 * data01$d20 
           + b40 * data01$d40 
           + b80 * data01$d80 
           - (b5 + b10 + b20 + b40 + b80) * data01$d160
           )

# update formals
formals(nll_basic3) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a5 = a5, a10 = a10, a20 = a20, a40 = a40, a80 = a80, 
  b2 = b2, b5 = b5, b10 = b10, b20 = b20, b40 = b40, b80 = b80,
  data = data01, 
  time = t, censor = censored, infected_treatment = g, 
  d1 = 'Gumbel', d2 = 'Weibull'
  )

m03_prep_function <- function(
  a1, b1, 
  a2, a5, a10, a20, a40, a80, 
  b2, b5, b10, b20, b40, b80){
    nll_basic3(
      a1, b1, 
      a2, a5, a10, a20, a40, a80, 
      b2, b5, b10, b20, b40, b80)
    }

m03 <- mle2(m03_prep_function, 
  start = list(
  a1 = 23, b1 = 4.6,
  a2 = 3, a5 = 0, a10 = 0, a20 = 0, a40 = 0, a80 = 0, 
  b2 = 0.2, b5 = 0, b10 = 0, b20 = 0, b40 = 0, b80 = 0)
  )

summary(m03)

### Model 4 
# this model estimates the full dose * food interaction
# and replaces the approximate approach used in the original text

data01 <- data_lorenz
head(data01)

# create columns with dummy variables

data01$d5 <- ifelse(data01$Infectious.dose == 5000, 1, 0)
data01$d10 <- ifelse(data01$Infectious.dose == 10000, 1, 0)
data01$d20 <- ifelse(data01$Infectious.dose == 20000, 1, 0)
data01$d40 <- ifelse(data01$Infectious.dose == 40000, 1, 0)
data01$d80 <- ifelse(data01$Infectious.dose == 80000, 1, 0)
data01$d160 <- ifelse(data01$Infectious.dose == 160000, 1, 0)

data01$af50 <- ifelse(data01$Food == 50, 1, 0)
data01$af100 <- ifelse(data01$Food == 100, 1, 0)

data01$f50d5 <- ifelse(data01$Infectious.dose == 5000 & data01$Food == 50, 1, 0)
data01$f50d10 <- ifelse(data01$Infectious.dose == 10000 & data01$Food == 50, 1, 0)
data01$f50d20 <- ifelse(data01$Infectious.dose == 20000 & data01$Food == 50, 1, 0)
data01$f50d40 <- ifelse(data01$Infectious.dose == 40000 & data01$Food == 50, 1, 0)
data01$f50d80 <- ifelse(data01$Infectious.dose == 80000 & data01$Food == 50, 1, 0)
data01$f50d160 <- ifelse(data01$Infectious.dose == 160000 & data01$Food == 50, 1, 0)

data01$f100d5 <- ifelse(data01$Infectious.dose == 5000 & data01$Food == 100, 1, 0)
data01$f100d10 <- ifelse(data01$Infectious.dose == 10000 & data01$Food == 100, 1, 0)
data01$f100d20 <- ifelse(data01$Infectious.dose == 20000 & data01$Food == 100, 1, 0)
data01$f100d40 <- ifelse(data01$Infectious.dose == 40000 & data01$Food == 100, 1, 0)
data01$f100d80 <- ifelse(data01$Infectious.dose == 80000 & data01$Food == 100, 1, 0)
data01$f100d160 <- ifelse(data01$Infectious.dose == 160000 & data01$Food == 100, 1, 0)

head(data01)

nll_basic4 <- nll_basic

# make 'pfa2' and 'pfb2' functions of food-by-dose interaction

body(nll_basic4)[[4]] <- substitute(
pfa2 <- a2 + a5 * data01$d5 
           + a10 * data01$d10 
           + a20 * data01$d20 
           + a40 * data01$d40 
           + a80 * data01$d80 
           - (a5 + a10 + a20 + a40 + a80) * data01$d160
           + af * data01$af50
           - af * data01$af100
           + afd5 * data01$f50d5
           + afd10 * data01$f50d10
           + afd20 * data01$f50d20
           + afd40 * data01$f50d40
           + afd80 * data01$f50d80
           - (afd5 + afd10 + afd20 + afd40 + afd80) * data01$f50d160
           - afd5 * data01$f100d5
           - afd10 * data01$f100d10
           - afd20 * data01$f100d20
           - afd40 * data01$f100d40
           - afd80 * data01$f100d80
           + (afd5 + afd10 + afd20 + afd40 + afd80) * data01$f100d160
           )

body(nll_basic4)[[5]] <- substitute(
pfb2 <- b2 + b5 * data01$d5 
           + b10 * data01$d10 
           + b20 * data01$d20 
           + b40 * data01$d40 
           + b80 * data01$d80 
           - (b5 + b10 + b20 + b40 + b80) * data01$d160
           + bf * data01$af50
           - bf * data01$af100
           + bfd5 * data01$f50d5
           + bfd10 * data01$f50d10
           + bfd20 * data01$f50d20
           + bfd40 * data01$f50d40
           + bfd80 * data01$f50d80
           - (bfd5 + bfd10 + bfd20 + bfd40 + bfd80) * data01$f50d160
           - bfd5 * data01$f100d5
           - bfd10 * data01$f100d10
           - bfd20 * data01$f100d20
           - bfd40 * data01$f100d40
           - bfd80 * data01$f100d80
           + (bfd5 + bfd10 + bfd20 + bfd40 + bfd80) * data01$f100d160
           )

formals(nll_basic4) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a5 = a5, a10 = a10, a20 = a20, a40 = a40, a80 = a80,
  af = af, afd5 = afd5, afd10 = afd10, afd20 = afd20, afd40 = afd40, afd80 = afd80,
  b2 = b2, b5 = b5, b10 = b10, b20 = b20, b40 = b40, b80 = b80,
  bf = bf, bfd5 = bfd5, bfd10 = bfd10, bfd20 = bfd20, bfd40 = bfd40, bfd80 = bfd80,
  data = data01, 
  time = t, 
  censor = censored,
  infected_treatment = g, 
  d1 = 'Gumbel', d2 = 'Weibull') 

m04_prep_function <- function(
  a1 = a1, b1 = b1, 
  a2 = a2, a5 = a5, a10 = a10, a20 = a20, a40 = a40, a80 = a80,
  af = af, afd5 = afd5, afd10 = afd10, afd20 = afd20, afd40 = afd40, afd80 = afd80,
  b2 = b2, b5 = b5, b10 = b10, b20 = b20, b40 = b40, b80 = b80,
  bf = bf, bfd5 = bfd5, bfd10 = bfd10, bfd20 = bfd20, bfd40 = bfd40, bfd80 = bfd80
  ){nll_basic4(
  a1 = a1, b1 = b1, 
  a2 = a2, a5 = a5, a10 = a10, a20 = a20, a40 = a40, a80 = a80,
  af = af, afd5 = afd5, afd10 = afd10, afd20 = afd20, afd40 = afd40, afd80 = afd80,
  b2 = b2, b5 = b5, b10 = b10, b20 = b20, b40 = b40, b80 = b80,
  bf = bf, bfd5 = bfd5, bfd10 = bfd10, bfd20 = bfd20, bfd40 = bfd40, bfd80 = bfd80
  )}

m04 <- mle2(m04_prep_function,
  start = list(
    a1 = 23, b1 = 4.6, 
    a2 = 3, a5 = 0.18, a10 = 0.03, a20 = 0.04, a40 = -0.05, a80 = -0.08,
    af = 0, afd5 = 0, afd10 = 0, afd20 = 0, afd40 = 0, afd80 = 0,
    b2 = 0.2, b5 = -0.03, b10 = 0.09, b20 = -0.01, b40 = 0.01, b80 = -0.02,
    bf = 0, bfd5 = 0, bfd10 = 0, bfd20 = 0, bfd40 = 0, bfd80 = 0)
  )

coef(m04)

### Model 5
# recall 'nll_basic3' from Model 3 above
# return 'pfb2' to original form 

body(nll_basic3) 

body(nll_basic3)[[5]] <- substitute(pfb2 <- b2)

formals(nll_basic3) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a5 = a5, a10 = a10, a20 = a20, a40 = a40, a80 = a80, 
  b2 = b2,
  data = data01, 
  time = t, censor = censored, infected_treatment = g, 
  d1 = 'Gumbel', d2 = 'Weibull'
  )

m05_prep_function <- function(a1, b1, a2, a5, a10, a20, a40, a80, b2){
    nll_basic3(a1, b1, a2, a5, a10, a20, a40, a80, b2)
    }

m05 <- mle2(m05_prep_function,
  start = list(
    a1 = 23, b1 = 4.6, 
    a2 = 3, a5 = 0.18, a10 = 0.03, a20 = 0.04, a40 = -0.05, a80 = -0.08, 
    b2 = 0.2)
  )

summary(m05)

### Model 6
# make 'pfa2' a linear function of log(Infectious dose)

nll_basic6 <- nll_basic 

body(nll_basic6)[[4]] <- substitute(pfa2 <- a2i + a2ii * log(data01$Infectious.dose))

formals(nll_basic6) <- alist(
  a1 = a1, b1 = b1, a2i = a2i, a2ii = a2ii, b2 = b2,
  data = data01, time = t, censor = censored, infected_treatment = g, 
  d1 = 'Gumbel', d2 = 'Weibull')

m06_prep_function <- function(a1, b1, a2i, a2ii, b2){
  nll_basic6(a1, b1, a2i, a2ii, b2)
  }

m06 <- mle2(m06_prep_function,
  start = list(a1 = 23, b1 = 4.6, a2i = 4, a2ii = -0.1, b2 = 0.2)
  )

summary(m06)

AICc(m01, m02, m03, m04, m05, m06, nobs = 256)


### Model 6 with different distributions 

m06b_prep_function <- function(a1, b1, a2i, a2ii, b2){
  nll_basic6(a1, b1, a2i, a2ii, b2, d1 = 'Weibull', d2 = 'Gumbel')
  }

m06b <- mle2(m06b_prep_function,
  start = list(a1 = 2, b1 = 0.5, a2i = 23, a2ii = -0.1, b2 = 4)
  )

m06c_prep_function <- function(a1, b1, a2i, a2ii, b2){
  nll_basic6(a1, b1, a2i, a2ii, b2, d1 = 'Weibull', d2 = 'Weibull')
  }

m06c <- mle2(m06c_prep_function,
  start = list(a1 = 2, b1 = 0.5, a2i = 3.8, a2ii = -0.1, b2 = 0.2)
  )

m06d_prep_function <- function(a1, b1, a2i, a2ii, b2){
  nll_basic6(a1, b1, a2i, a2ii, b2, d1 = 'Gumbel', d2 = 'Gumbel')
  }

m06d <- mle2(m06d_prep_function,
  start = list(a1 = 23, b1 = 5, a2i = 23, a2ii = -0.1, b2 = 5)
  )

AICc(m06, m06b, m06c, m06d, nobs = 256)


## ----warning = FALSE----------------------------------------------------------

 # data01 <- data_blanford_bl5
 # subset data for 'block = 5', treatments 'cont', 'Ma06', 'Ma07', 'Ma08'

data01 <- data_blanford 

data01 <- subset(data01,
    (data01$block == 5) & (
      (data01$treatment == 'cont') |
      (data01$treatment == 'Ma06') |
      (data01$treatment == 'Ma07') |
      (data01$treatment == 'Ma08') ) &
    (data01$day > 0)
    )

# create column 'g' as index of infected treatment
data01$g <- data01$inf

head(data01)

nll_basic2 <- nll_basic

# make 'pfa2' and 'pfb2' functions of fungal treatment
# NB to avoid problems with log(0), set final ifelse 'false' values to 'exp(0)'

body(nll_basic2)[[4]] <- substitute(pfa2 <- 
  ifelse(((data01$g == 1) & (data01$treatment == 'Ma06')), a2, 
  ifelse(((data01$g == 1) & (data01$treatment == 'Ma07')), a3,
  ifelse(((data01$g == 1) & (data01$treatment == 'Ma08')), a4,
  exp(0)
  ))))

body(nll_basic2)[[5]] <- substitute(pfb2 <- 
  ifelse(((data01$g == 1) & (data01$treatment == 'Ma06')), b2, 
  ifelse(((data01$g == 1) & (data01$treatment == 'Ma07')), b3,
  ifelse(((data01$g == 1) & (data01$treatment == 'Ma08')), b4,
  exp(0)
  ))))

formals(nll_basic2) <- alist(
  a1 = a1, b1 = b1,
  a2 = a2, b2 = b2,
  a3 = a3, b3 = b3,
  a4 = a4, b4 = b4, 
  data = data01, 
  time = t, censor = censor, infected_treatment = g,
  d1 = 'Weibull', d2 = 'Fréchet')

m01_prep_function <- function(a1, b1, a2, b2, a3, b3, a4, b4){
  nll_basic2(a1, b1, a2, b2, a3, b3, a4, b4)
  }

m01 <- mle2(m01_prep_function, 
  start = list(a1 = 2, b1 = 1, a2 = 2, b2 = 1, a3 = 2, b3 = 1, a4 = 2, b4 = 1)
  )

summary(m01)
confint(m01)


## ----warning = FALSE----------------------------------------------------------

data01 <- data_parker

head(data01, 3)

### Model 1 

nll_basic2 <- nll_basic

body(nll_basic2)[[4]] <- substitute(pfa2 <- 
  a2 + ifelse(data01$dose == 1, a2i, 
       ifelse(data01$dose == 2, a2ii, 
       ifelse(data01$dose == 3, -(a2i + a2ii), 
       exp(0)
       )))
  )

body(nll_basic2)[[5]] <- substitute(pfb2 <- 
  b2 + ifelse(data01$dose == 1, b2i, 
       ifelse(data01$dose == 2, b2ii, 
       ifelse(data01$dose == 3, -(b2i + b2ii), 
       exp(0)
       )))
  )

formals(nll_basic2) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a2i = a2i, a2ii = a2ii, 
  b2 = b2, b2i = b2i, b2ii = b2ii,
  data = data01,
  time = t, censor = censored, infected_treatment = g,
  d1 = 'Frechet', d2 = 'Frechet')

m01_prep_function <- function(
  a1, b1, a2, a2i, a2ii, b2, b2i, b2ii){
    nll_basic2(a1, b1, a2, a2i, a2ii, b2, b2i, b2ii)
  }

m01 <- mle2(m01_prep_function,
  start = list(
    a1 = 2, b1 = 1,
    a2 = 2, a2i = 0, a2ii = 0, 
    b2 = 1, b2i = 0, b2ii = 0)
  )

summary(m01)

### Model 2
# make 'pfa2' and 'pfb2' linear functions of log(dose)

body(nll_basic2)[[4]] <- substitute(pfa2 <- 
  a2 + ifelse(data01$dose == 1, a2i, 
       ifelse(data01$dose == 2, 0, 
       ifelse(data01$dose == 3, -a2i, 
       exp(0)
       )))
  )

body(nll_basic2)[[5]] <- substitute(pfb2 <- 
  b2 + ifelse(data01$dose == 1, b2i, 
       ifelse(data01$dose == 2, 0, 
       ifelse(data01$dose == 3, -b2i, 
       exp(0)
       )))
  ) 

formals(nll_basic2) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a2i = a2i, 
  b2 = b2, b2i = b2i,
  data = data01,
  time = t, censor = censored, infected_treatment = g,
  d1 = 'Frechet', d2 = 'Frechet')

m02_prep_function <- function(
  a1, b1, a2, a2i, b2, b2i){
    nll_basic2(a1, b1, a2, a2i, b2, b2i)
  }

m02 <- mle2(m02_prep_function,
  start = list(
    a1 = 2, b1 = 1, a2 = 2, a2i = 0, b2 = 1, b2i = 0)
    )

summary(m02)

### Model 3

nll_basic3 <- nll_basic

body(nll_basic3)[[4]] <- substitute(pfa2 <- 
  a2 + ifelse(data01$g == 1, 
       ifelse(data01$Sporulation == 1, a2sp, -a2sp), 
       exp(0))
  )

body(nll_basic3)[[5]] <- substitute(pfb2 <- 
  b2 + ifelse(data01$g == 1, 
       ifelse(data01$Sporulation == 1, b2sp, -b2sp), 
       exp(0))
  )

formals(nll_basic3) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a2sp = a2sp, 
  b2 = b2, b2sp = b2sp,
  data = data01,
  time = t, censor = censored, infected_treatment = g,
  d1 = 'Fréchet', d2 = 'Weibull')

m03_prep_function <- function(a1, b1, a2, a2sp, b2, b2sp){
  nll_basic3(a1, b1, a2, a2sp, b2, b2sp)
  }

m03 <- mle2(m03_prep_function,
  start = list(
    a1 = 2.35, b1 = 0.66, 
    a2 = 2, a2sp = 0,
    b2 = 0.5, b2sp = 0)
  )

summary(m03) 

### Model 4

nll_basic4 <- nll_basic

data01 <- data_parker

body(nll_basic4)[[4]] <- substitute(pfa2 <- 
  a2 + ifelse(data01$g == 1, 
       ifelse(data01$Sporulation == 1, 
       ifelse(data01$dose == 1, a2sp + a2d1, 
       ifelse(data01$dose == 3, a2sp - a2d1, a2sp)), 
       exp(0)), 
       exp(0))
  )

body(nll_basic4)[[5]] <- substitute(pfb2 <- 
  b2 + ifelse(data01$g == 1, 
       ifelse(data01$Sporulation == 1, 
       ifelse(data01$dose == 1, b2sp + b2d1, 
       ifelse(data01$dose == 3, b2sp - b2d1, b2sp)), 
       exp(0)), 
       exp(0))
  )

formals(nll_basic4) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a2sp = a2sp, a2d1 = a2d1,
  b2 = b2, b2sp = b2sp, b2d1 = b2d1,
  data = data01,
  time = t, censor = censored, infected_treatment = g,
  d1 = 'Fréchet', d2 = 'Weibull')

m04_prep_function <- function(a1, b1, a2, a2sp, a2d1, b2, b2sp, b2d1){
    nll_basic4(a1, b1, a2, a2sp, a2d1, b2, b2sp, b2d1)
  } 

m04 <- mle2(m04_prep_function,
  start = list(
    a1 = 2.35, b1 = 0.66,
    a2 = 2, a2sp = 0, a2d1 = 0,
    b2 = 0.5, b2sp = 0, b2d1 = 0)
  )

summary(m04)

### Model 5
 
# NB here analyse sporulating vs unsporulating hosts
  # i.e. pool uninfected + sporulation = 0 hosts together
  # so for index of 'infected_treatment' use 'Sporulation'

nll_basic5 <- nll_basic

body(nll_basic5)[[4]] <- substitute(pfa2 <- 
  ifelse(data01$dose == 1, a2 + a2d1, 
  ifelse(data01$dose == 3, a2 - a2d1, 
  a2))
  )

formals(nll_basic5) <- alist(
  a1 = a1, b1 = b1, 
  a2 = a2, a2d1 = a2d1, b2 = b2,
  data = data01,
  time = t, 
  censor = censored, 
  infected_treatment = Sporulation,
  d1 = 'Fréchet', d2 = 'Weibull')


m05_prep_function <- function(a1, b1, a2, a2d1, b2){
  nll_basic5(a1, b1, a2, a2d1, b2)
  }

m05 <- mle2(m05_prep_function, 
  start = list(a1 = 2, b1 = 1, a2 = 2, a2d1 = 0.1, b2 = 0.5)
  )

summary(m05)

AICc(m01, m02, m03, m04, m05, nobs = 328)


## ----warning = FALSE----------------------------------------------------------

data01 <- data_parker

head(data01)

# Infection status known 

# Here a host's infection status is defined by whether
# it had visible signs of sporulation at the time of its death
# or right-censoring 
# non-sporulating hosts are assumed to only experience same background
# mortality as uninfected hosts

# this is equivalent to taking infection_treatment =  Sporulation


m01_prep_function <- function(a1, b1, a2, b2){
  nll_basic(a1, b1, a2, b2,
    data = data01,
    time = t,
    censor = censored, 
    infected_treatment = Sporulation,
    d1 = 'Fréchet', d2 = 'Weibull')
  }

m01 <- mle2(m01_prep_function, 
  start = list(a1 = 2.56, b1 = 0.72, a2 = 2, b2 = 0.5)
  )

summary(m01)


# Infection status unknown 

# the infection status of individual hosts is not always known
# the observed survival data may suggest an infected treatment 
# harboured exposed-infected and exposed-uninfected hosts
# nll_exposed_infected estimates the proportion of hosts in
# an infected treatment experiencing increased rates of mortality 
# due to infection (p1), and
# the proportion experiencing only background mortality (1 - p1)

m02_prep_function <- function(a1, b1, a2, b2, p1){
  nll_exposed_infected(a1, b1, a2, b2, p1,
    data = data01, 
    time = t,
    censor = censored,
    infected_treatment = g, 
    d1 = 'Frechet', d2 = 'Weibull')
  }

m02 <- mle2(m02_prep_function, 
  start = list(a1 = 2, b1 = 1, a2 = 2, b2 = 0.5, p1 = 0.5)
  ) 

summary(m02)

# in the Parker et al data, the proportion of sporulating hosts in the
# infected treatments were, 119 / (119 + 127) = 0.484

aggregate(data01, by = list(data01$g, data01$Sporulation), length)


# extend the above to the proportion sporulating within dose treatments

nll_exposed_infected2 <- nll_exposed_infected

body(nll_exposed_infected2)[[6]] <- substitute(pfp1 <- 
  p1 + ifelse(data01$g == 1 & data01$dose == 1, -p1d,
       ifelse(data01$g == 1 & data01$dose == 3, + p1d, 
       ifelse(data01$g == 1 & data01$dose == 2, 0,
       exp(0)
       )))
  )

formals(nll_exposed_infected2) <- alist(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2,
  p1 = p1, p1d = p1d,
  data = data01,
  time = t,
  censor = censored,
  infected_treatment = g,
  d1 = 'Frechet', d2 = 'Weibull'
  )

m03_prep_function <- function(a1, b1, a2, b2, p1, p1d){
  nll_exposed_infected2(a1, b1, a2, b2, p1, p1d)
  }

m03 <- mle2(m03_prep_function,
  start = list(a1 = 2.2, b1 = 0.53, a2 = 1.88, b2 = 0.17, 
    p1 = 0.48, p1d = 0)
  )

summary(m03)

# estimated proportions infected;
  # dose 1 = 0.51 - 0.24 = 0.27
  # dose 2 = 0.51 + 0    = 0.51
  # dose 3 = 0.51 + 0.24 = 0.75

# observed proportions sporulating;

  aggregate(data01, by = list(data01$g, data01$Sporulation, data01$dose), length)

# dose 1 = 25 / (25 + 56) = 0.31
# dose 2 = 40 / (40 + 42) = 0.49
# dose 3 = 54 / (54 + 29) = 0.65
 
AICc (m01, m02, m03, nobs = 328)


## ----warning = FALSE----------------------------------------------------------

data01 <- data_lorenz

m01_prep_function <- function(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2){
      nll_basic(
        a1 = a1, b1 = b1, a2 = a2, b2 = b2,
        data = data01,
        time = t,
        censor = censored,
        infected_treatment = g,
        d1 = "Gumbel",
        d2 = "Frechet"
        )}

m01 <- mle2(m01_prep_function,
      start = list(a1 = 20, b1 = 5, a2 = 3, b2 = 0.5)
      )

summary(m01)

m02_prep_function <- function(
  a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta){
      nll_frailty(
        a1 = a1, b1 = b1, a2 = a2, b2 = b2, theta = theta,
        data = data_lorenz,
        time = t,
        censor = censored,
        infected_treatment = g,
        d1 = "Gumbel",
        d2 = "Weibull",
        d3 = "Gamma"
        )}

m02 <- mle2(m02_prep_function,
      start = list(a1 = 20, b1 = 5, a2 = 3, b2 = 0.1, theta = 2)
      )

summary(m02)

AICc(m01, m02, nobs = 256) 

## ----warning = FALSE----------------------------------------------------------

data01 <- data_lorenz

m01_prep_function <- function(a1, b1, a2, b2, theta){
  nll_frailty_shared(a1, b1, a2, b2, theta,
    data = data01,
    time = t,
    censor = censored,
    infected_treatment = g,
    d1 = "Gumbel", d2 = "Gumbel")
  }


m01 <- mle2(m01_prep_function,
  start = list(a1 = 23, b1 = 5, a2 = 10, b2 = 1, theta = 1),
  method = "Nelder-Mead",
  control = list(maxit = 5000)
  )

summary(m01)


m02_prep_function <- function(a1, b1, a2, b2, theta01, theta02, rho){
  nll_frailty_correlated(a1, b1, a2, b2, theta01, theta02, rho,
    data = data01,
    time = t,
    censor = censored,
    infected_treatment = g,
    d1 = "Gumbel",
    d2 = "Gumbel")
  }

m02 <- mle2(m02_prep_function,
  start = list(
    a1 = 20, b1 = 5, a2 = 20, b2 = 4, theta01 = 1, theta02 = 1, rho = 1),
    method = "L-BFGS-B",
    lower = list(
      a1 = 1e-6, b1 = 1e-6, a2 = 1e-6, b2 = 1e-6,
      theta01 = 1e-6, theta02 = 1e-6, rho = 1e-6)
    )

summary(m02)

# NB no standard errors estimated and estimate of 'theta01' at lower boundary
# rerun model with theta01 set at lower boundary
    
m02b <- mle2(m02_prep_function,
  start = list(
    a1 = 20, b1 = 5, a2 = 20, b2 = 4,
    theta01 = 1, theta02 = 1, rho = 1),
    fixed = list(theta01 = 1e-6),
    method = "L-BFGS-B",
    lower = list(
      a1 = 1e-6, b1 = 1e-6, a2 = 1e-6, b2 = 1e-6, 
      theta02 = 1e-6, rho = 1e-6)
    )

summary(m02b)

# NB standard error of 'rho' crosses zero (0)
# rerun model with rho set to lower limit

m02c <- mle2(m02_prep_function,
  start = list(
    a1 = 20, b1 = 5, a2 = 20, b2 = 4,
    theta01 = 1, theta02 = 1, rho = 1),
    fixed = list(theta01 = 1e-6, rho = 1e-6),
    method = "L-BFGS-B",
    lower = list(
      a1 = 1e-6, b1 = 1e-6, a2 = 1e-6, b2 = 1e-6, 
      theta02 = 1e-6)
    )

summary(m02c)

# result of m02c corresponds with estimates from the 'nll_frailty' model,
# where it is assumed there is no unobserved variation in the rate of background mortality 
# and where the gamma distribution describes the unobserved variation in virulence

m03_prep_function <- function(a1, b1, a2, b2, theta){
      nll_frailty(a1, b1, a2, b2, theta,
        data = data01, time = t,
        censor = censored, infected_treatment = g,
        d1 = "Gumbel", d2 = "Gumbel", d3 = "Gamma")
      }

m03 <- mle2(m03_prep_function,
  start = list(a1 = 20, b1 = 4, a2 = 20, b2 = 4, theta = 1)
  )

summary(m03)

AICc(m02c, m03, nobs = 256)

coef(m02c)
coef(m03)


