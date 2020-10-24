## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(anovir)

## ----include = FALSE----------------------------------------------------------

sy_time <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)

sy_Suninf <- c(1,0.988,0.988,0.969,0.951,0.941,0.934,0.929,0.92,0.913,0.906,0.903,0.896,0.889,0.88,0.875,0.875,0.868,0.858,0.846,0.842,0.837,0.837,0.831)

sy_Sobsinf <- c(1,0.979,0.966,0.942,0.915,0.902,0.887,0.875,0.869,0.859,0.841,0.81,0.788,0.747,0.732,0.696,0.68,0.652,0.616,0.6,0.584,0.574,0.561,0.545)

sy_surv <- data.frame(sy_time, sy_Suninf, sy_Sobsinf)


## ----echo = FALSE, fig.width = 4, fig.height = 4, fig.align = 'center'--------

   plot(sy_surv[, 1], sy_surv[, 2],
     type = 's', col = 'black', xlab = 'time (days)', ylab = 'Survival',
     xlim = c(0, 23), ylim = c(0, 1),
     main = 'Cumulative survival',
     axes = FALSE
     )
     lines(sy_surv[, 3], type = 's', col = 'red')
     axis(side = 1, seq(0, 28, by = 7))
     axis(side = 2, seq(0, 1.1, by = 0.2))
     legend("bottomleft", c("Uninfected", "Infected"), lty = c(1, 1), col = c(1, 2))


## ----include = FALSE----------------------------------------------------------

  sy_times <- c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,23,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,7,21,22,23,7,14,21,22,23)
  sy_censor <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
  sy_inf <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1)
  sy_fq <- c(5,8,8,4,3,2,4,3,3,1,3,3,4,2,3,4,5,2,2,1,7,4,8,9,4,5,4,2,3,6,10,7,13,5,11,5,9,11,5,5,3,4,2,3,4,205,143,6,7,4,103,66)

  data_sy <- data.frame(sy_times, sy_inf, sy_censor, sy_fq)
  
  data_sy$fq <- data_sy$sy_fq

## ----include = TRUE-----------------------------------------------------------
    head(data_sy)

## ----eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE---------------

    m01_prep_function <- function(a1, b1, a2, b2){
      nll_basic(
        a1, b1, a2, b2,
        data = data_sy,
        time = sy_times,
        censor = sy_censor,
        infected_treatment = sy_inf,
        d1 = 'Weibull', d2 = 'Weibull')
        }

    m01 <- mle2(m01_prep_function,
             start = list(a1 = 3, b1 = 0.5, a2 = 3, b2 = 0.5)
             )
    
    summary(m01)
    

## ----eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE---------------
    m02 <- mle2(m01_prep_function,
             start = list(a1 = 4.8, a2 = 3.6, b2 = 0.6),
             fixed = list(b1 = 1.0)
             )
    summary(m02)

  # compare models by AICc
     AICc(m01, m02, nobs = 753)

