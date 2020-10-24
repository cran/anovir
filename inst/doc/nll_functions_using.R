## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(anovir)

## ----include = FALSE----------------------------------------------------------
data01 <- subset(data_blanford,
  (data_blanford$block == 3) & 
  ((data_blanford$treatment == 'cont') | (data_blanford$treatment == 'Bb06')) &
  (data_blanford$day > 0)
  )

## -----------------------------------------------------------------------------
head(data01, 6)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#      m01_prep_function <- function(a1, b1, a2, b2){
#        nll_basic(
#          a1, b1, a2, b2,
#          data = data01,
#          time = day,
#          censor = censor,
#          infected_treatment = inf,
#          d1 = 'Weibull', d2 = 'Weibull')
#          }

## ----eval = FALSE, echo = TRUE------------------------------------------------
#      m01 <- mle2(m01_prep_function,
#               start = list(a1 = 2, b1 = 0.5, a2 = 2, b2 = 0.5)
#               )

## ----eval = TRUE, echo = FALSE, collapse = TRUE, error = FALSE, message = FALSE, warning = FALSE----

    m01_prep_function <- function(a1, b1, a2, b2){
      nll_basic(
        a1, b1, a2, b2,
        data = data01,
        time = day,
        censor = censor,
        infected_treatment = inf,
        d1 = 'Weibull',
        d2 = 'Weibull')
        }

    m01 <- mle2(m01_prep_function,
             start = list(a1 = 2, b1 = 0.5, a2 = 2, b2 = 0.5)
             )

    summary(m01)


## ----fig.width = 4, fig.height = 4, fig.align = 'center'----------------------

    coef(m01)
    
    vcov(m01)

    a2 <- coef(m01)[[3]]
    b2 <- coef(m01)[[4]]
    var_a2 <- vcov(m01)[3,3]
    var_b2 <- vcov(m01)[4,4]
    cov_a2b2 <- vcov(m01)[3,4]
    
    ci_matrix01 <- conf_ints_virulence(
      a2 = a2, b2 = b2,
      var_a2 = var_a2, var_b2 = var_b2, cov_a2b2 = cov_a2b2, 
      d2 = "Weibull", tmax = 14)

    tail(ci_matrix01)

    plot(ci_matrix01[, 't'], ci_matrix01[, 'h2'], 
         type = 'l', col = 'red', xlab = 'time', ylab = 'virulence (± 95% ci)'
         )
      lines(ci_matrix01[, 'lower_ci'], col = 'grey')
      lines(ci_matrix01[, 'upper_ci'], col = 'grey')

