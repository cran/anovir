## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(anovir)

## ----include = FALSE----------------------------------------------------------
subset_data_bl3 <- subset(data_blanford,
  (data_blanford$block == 3) & 
  ((data_blanford$treatment == 'cont') | (data_blanford$treatment == 'Bb06')) &
  (data_blanford$day > 0)
  )

## ----include = FALSE, message = FALSE, warning = FALSE------------------------
  time <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
  Suninf <- c(0.976,0.96,0.952,0.944,0.927,0.911,0.879,0.831,0.782,0.742,0.726,0.685,0.613,0.508)
  Sobsinf <- c(0.988,0.969,0.963,0.939,0.902,0.853,0.779,0.687,0.601,0.509,0.399,0.337,0.172,0.11)
  data01 <- data.frame(time, Suninf, Sobsinf)
  data01$log_time <- round(log(data01$time), 3)
  data01$clogSuninf <- round(log(-log(data01$Suninf)), 3)
  data01$clogSobsinf <- round(log(-log(data01$Sobsinf)), 3)
  data01$Srel <- data01$Sobsinf / data01$Suninf
  data01$clogSrel <- log(-log(data01$Srel))

## ---- echo = FALSE, fig.height = 3.25, fig.width = 3.25, fig.show = 'hold'----
    plot(data01$time, data01$Suninf, 
       type = 's', col = 'black', lty = 'solid',
       xlim = c(0, 15), ylim = c(0, 1.1),
       xlab = 'time (days)', ylab = 'Cumulative survival',
       main = 'Observed survival',
       axes = FALSE
       )
  lines(data01$time, data01$Suninf, 
        type = 'p', col = 'black', cex = 0.5
        )
  lines(data01$time, data01$Sobsinf, 
        type = 's', col = 'blue', lty = 'solid'
        )
  lines(data01$time, data01$Sobsinf, 
        type = 'p', col = 'blue', cex = 0.5
        )
    axis(side = 1, seq(0, 15, by = 7))
    axis(side = 2, seq(0, 1.1, by = 0.25))
#    legend("bottomleft", c("Uninfected", "Infected"), lty = c(1, 1), col = c('black', 'blue'))
  
  plot(data01$log_time, data01$clogSuninf, 
       type = 'l', col = 'black', lty = 'solid',
       xlim = c(0, 3), ylim = c(-5, 1),
       xlab = 'log(time)', ylab = 'log(-log S[t])',
       main = 'Transformed survival',
       axes = FALSE
       )
  lines(data01$log_time, data01$clogSuninf, 
        type = 'p', col = 'black', cex = 0.5
        )
  lines(data01$log_time, data01$clogSobsinf, 
        type = 'l', col = 'blue', lty = 'solid'
        )
  lines(data01$log_time, data01$clogSobsinf, 
        type = 'p', col = 'blue', cex = 0.5
        )
    axis(side = 1, seq(0, 3, by = 1))
    axis(side = 2, seq(-5, 1, by = 2))
#    legend("bottomleft", c("Uninfected", "Infected"), lty = c(1, 1), col = c('black', 'blue'))

## ---- include = TRUE----------------------------------------------------------
lr_clogSuninf <- lm(data01$clogSuninf ~ data01$log_time, data = data01)
coef(lr_clogSuninf)
df.residual(lr_clogSuninf)

## ---- include = TRUE----------------------------------------------------------
lr_clogSrel <- lm(data01$clogSrel ~ data01$log_time, data = data01)
coef(lr_clogSrel)
df.residual(lr_clogSrel)

## ---- message = FALSE, warning = FALSE----------------------------------------
     nlr01 <- nls(Suninf ~ exp(-exp((log(time) - a1) / b1)), 
                data = data01,
                start = list(a1 = 2, b1 = 1)
                )
    coef(nlr01)
    df.residual(nlr01)

## ---- message = FALSE, warning = FALSE----------------------------------------
      nlr02 <- nls(Sobsinf ~ exp(-exp((log(time) - 2.888) / 0.474)) * 
                           exp(-exp((log(time) - a2) / b2)),
                           data = data01,
                           start = list(a2 = 2, b2 = 1)
                           )
      coef(nlr02)
      df.residual(nlr02)


## ----include = TRUE-----------------------------------------------------------
  head(subset_data_bl3, 3)

## ---- include = TRUE, message = FALSE, warning = FALSE------------------------
  # step #1: prep-function for maximum likelihood estimation
     mle01_prep_function <- function(a1, b1, a2, b2){
       nll_basic(
         a1, b1, a2, b2,
         data = subset_data_bl3,
         time = day,
         censor = censor,
         infected_treatment = inf,
         d1 = 'Weibull', d2 = 'Weibull'
       )}

  # step #2: start values from linear regression of transformed data
    mle01 <- mle2(mle01_prep_function,
               start = list(a1 = 3.339, b1 = 0.790, a2 = 2.515, b2 = 0.239)
               )
   coef(mle01)


## ---- include = TRUE----------------------------------------------------------
   # linear regression estimates
    LR <- mle2(mle01_prep_function,
                 start = list(a1 = 3.339, b1 = 0.790, a2 = 2.515, b2 = 0.239),
                 eval.only = TRUE
                 )

   # non-linear regression estimates   
   NLR <- mle2(mle01_prep_function,
                 start = list(a1 = 2.888, b1 = 0.474, a2 = 2.542, b2 = 0.266),
                 eval.only = TRUE
                 )
   
   # maximum likelihood estimates
   MLE <- mle2(mle01_prep_function,
                 start = list(a1 = 2.845, b1 = 0.483, a2 = 2.581, b2 = 0.183),
                 eval.only = TRUE
                 )
   
   # compare models by Akaike's Information Criteria, adjusted for small sample sizes
   AICc(LR, NLR, MLE, nobs = 287)
   

## ---- include = TRUE----------------------------------------------------------
  coef(mle01) ; vcov(mle01)

  ci_matrix01 <- conf_ints_virulence(
    a2 = 2.5807774,
    b2 = 0.1831184,
    var_a2 = 0.0008196247,
    var_b2 = 0.0010005684,
    cov_a2b2 = -0.0003118982, 
    d2 = "Weibull", tmax = 14)


## ---- echo = FALSE, fig.height = 3.25, fig.width = 3.25, fig.show = 'hold'----
  a1 <- 2.845
  b1 <- 0.483
  a2 <- 2.581
  b2 <- 0.183

  data01$Suninf_est <- exp(-exp((log(time) - a1) / b1))
  data01$Sobsinf_est <- exp(-exp((log(time) - a1) / b1)) * exp(-exp((log(time) - a2) / b2))

  
  plot(data01$time, data01$Suninf, 
       type = 's', col = 'black', lty = 'dotted',
       xlim = c(0, 15), ylim = c(0, 1.1),
       xlab = 'time (days)', ylab = 'Cumulative survival',
       main = 'Estimated curves',
       axes = FALSE
       )
  lines(data01$time, data01$Suninf, 
        type = 'p', col = 'black', cex = 0.5
        )
  lines(data01$time, data01$Sobsinf, 
        type = 's', col = 'blue', lty = 'dotted'
        )
  lines(data01$time, data01$Sobsinf, 
        type = 'p', col = 'blue', cex = 0.5
        )
  lines(data01$Suninf_est, lty = 'solid', col = 'black', lwd = 2)
  lines(data01$Sobsinf_est, lty = 'solid', col = 'blue', lwd = 2)
    axis(side = 1, seq(0, 15, by = 7))
    axis(side = 2, seq(0, 1.1, by = 0.25))
 
  ci_matrix01 <- conf_ints_virulence(
    a2 = 2.5807774,
    b2 = 0.1831184,
    var_a2 = 0.0008196247,
    var_b2 = 0.0010005684,
    cov_a2b2 = -0.0003118982, 
    d2 = "Weibull", tmax = 14) 
    
  
    plot(ci_matrix01[, 't'], ci_matrix01[, 'h2'], 
         xlim = c(0, 14), ylim = c(0, 1.1),
         type = 'l', col = 'blue', lwd = 2,
         xlab = 'time (days)', ylab = 'virulence (± 95% ci)',
         main = 'Estimated virulence',
         axes = FALSE
         )
      lines(ci_matrix01[, 'lower_ci'], col = 'grey')
      lines(ci_matrix01[, 'upper_ci'], col = 'grey')
    axis(side = 1, seq(0, 14, by = 7))
    axis(side = 2, seq(0, 1, by = 0.2))


