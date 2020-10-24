## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(anovir)

## -----------------------------------------------------------------------------
utils::str(nll_basic)

## -----------------------------------------------------------------------------
head(body(nll_basic), 5)

## ----message = FALSE, warning = FALSE-----------------------------------------
    head(data_lorenz, 2)

    # step #1
       m01_prep_function <- function(a1 = a1, b1 = b1, a2 = a2, b2 = b2){
         nll_basic(a1 = a1, b1 = b1, a2 = a2, b2 = b2,
           data = data_lorenz, time = t, censor = censored,
           infected_treatment = g, d1 = 'Gumbel', d2 = 'Weibull')
       }
    # step #2
       m01 <- mle2(m01_prep_function,
                start = list(a1 = 23, b1 = 5, a2 = 3, b2 = 0.2)
                )
       coef(m01)

## ----message = FALSE, warning = FALSE-----------------------------------------
  # copy/rename 'nll_function' (not obligatory, but recommended)
    nll_basic2 <- nll_basic
  # find/check location of code to be replaced. NB double '[['
    body(nll_basic2)[[4]]
  # replace default code with new code for function
    body(nll_basic2)[[4]] <- substitute(pfa2 <- c1 + c2 * log(data$Infectious.dose))
  # check code
    head(body(nll_basic2), 5)
  # replace argument 'a2' with those for 'c1', 'c2'. NB use of 'alist'
    formals(nll_basic2) <- alist(a1 = a1, b1 = b1, c1 = c1, c2 = c2, b2 = b2,
                            data = data, time = time, censor = censor,
                            infected_treatment = infected_treatment, d1 = "", d2 = "")
  # new analysis: step #1
    m02_prep_function <- function(a1 = a1, b1 = b1, c1 = c1, c2 = c2, b2 = b2){
      nll_basic2(a1 = a1, b1 = b1, c1 = c1, c2 = c2, b2 = b2,
           data = data_lorenz, time = t, censor = censored,
           infected_treatment = g, d1 = 'Gumbel', d2 = 'Weibull')
         }
  # step #2
    m02 <- mle2(m02_prep_function,
                start = list(a1 = 23, b1 = 5, c1 = 4, c2 = -0.1, b2 = 0.2)
                )
    coef(m02)
    
  # compare results
    AICc(m01, m02, nobs = 256)
  # according to AICc m02 is better than m01  
    

## ----message = FALSE, warning = FALSE-----------------------------------------
  # copy/rename nll_function
    nll_basic3 <- nll_basic
    body(nll_basic3)[[4]] <- substitute(pfa2 <- c1 + c2 * log(data$Infectious.dose))

  # replace argument 'a2' with those for 'c1', 'c2', and assign column names  
    formals(nll_basic3) <- alist(a1 = a1, b1 = b1, c1 = c1, c2 = c2, b2 = b2,
                            data = data_lorenz, time = t, censor = censored,
                            infected_treatment = g, d1 = "Gumbel", d2 = "Weibull")
  # new analysis: step #1
    m03_prep_function <- function(a1 = a1, b1 = b1, c1 = c1, c2 = c2, b2 = b2){
      nll_basic3(a1 = a1, b1 = b1, c1 = c1, c2 = c2, b2 = b2)
      }
  # step #2
    m03 <- mle2(m03_prep_function,
                start = list(a1 = 23, b1 = 5, c1 = 4, c2 = -0.1, b2 = 0.2)
                )
    coef(m03)
    identical(coef(m02), coef(m03))


