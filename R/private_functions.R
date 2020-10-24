
# nb to pass check need to replace 'é' with '\u00E9'

# Private function, only for package
# Get zx as function of distribution type
# NB defaults to Weibull
P_get_zx <- function(t, ax, bx, dx){
  zx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    zx <- (log(t) - ax) / bx
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    zx <- (t - ax) / bx
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    zx <- (log(t) - ax) / bx
  } else {
    zx <- (log(t) - ax) / bx
  }
  return(zx)
}


# Private function, only for package
# Get fx as function of distribution type
P_get_fx <- function(t, zx, bx, dx){
  fx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    fx <- (1 / (bx * t)) * exp(zx - exp(zx))
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    fx <- (1 / bx) * exp(zx - exp(zx))
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    fx <- (1 / (bx * t)) * exp(-zx - exp(-zx))
  }  else {
    fx <- (1 / (bx * t)) * exp(zx - exp(zx))
  }
  return(fx)
}

# Private function, only for package
# Get hx as function of distribution type
P_get_hx <- function(t, zx, bx, dx){
  hx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    hx <- (1 / (bx * t)) * exp(zx)
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    hx <- (1 / bx) * exp(zx)
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    hx <- (1 / (bx * t)) * exp(-zx - exp(-zx)) / (1 - exp(-exp(-zx)))
  }  else {
    hx <- (1 / (bx * t)) * exp(zx)
  }
  return(hx)
}

# Private function, only for package
# Get Sx as function of distribution type
P_get_Sx <- function(t, zx, dx){
  Sx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    Sx <- exp(-exp(zx))
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    Sx <- exp(-exp(zx))
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    Sx <- 1 - exp(-exp(-zx))
  }  else {
    Sx <- exp(-exp(zx))
  }
  return(Sx)
}

# Private function, only for package
# Get Hx as function of distribution type
# to be used in estimating expected time of death (ETD) when Hx = 1
P_get_Hx <- function(t, zx, dx = ""){
  Hx <- 1
  if((dx == "Weibull") | (dx == "weibull")){
    Hx <- exp(zx)
  } else if ((dx == "Gumbel") | (dx == "gumbel")){
    Hx <- exp(zx)
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    Hx <- -log(1 - exp(-exp(-zx)))
  } else {
    Hx <- exp(zx)
  }
  return(Hx)
}


# Private function, only for package
# Get zu as function of distribution type
# NB defaults to Weibull
P_get_zu <- function(u, ax, bx, dx){
  zu <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    zu <- (log(u) - ax) / bx
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    zu <- (u - ax) / bx
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    zu <- (log(u) - ax) / bx
  } else {
    zu <- (log(u) - ax) / bx
  }
  return(zu)
}

# Private function, only for package
# Get Su as function of distribution type
P_get_Su <- function(u, zu, dx){
  Su <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    Su <- exp(-exp(zu))
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    Su <- exp(-exp(zu))
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    Su <- 1 - exp(-exp(-zu))
  }  else {
    Su <- exp(-exp(zu))
  }
  return(Su)
}

# Private function, only for package
# Get fu as function of distribution type
# used in recovery model
P_get_fu <- function(u, zu, bx, dx){
  fu <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    fu <- (1 / (bx * u)) * exp(zu - exp(zu))
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    fu <- (1 / bx) * exp(zu - exp(zu))
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    fu <- (1 / (bx * u)) * exp(-zu - exp(-zu))
  }  else {
    fu <- (1 / (bx * u)) * exp(zu - exp(zu))
  }
  return(fu)
}



# Private function, only for package
# Weibull functions for; zx, Sx, fx, hx, Hx in terms of zx
P_Weibull_functions_zx <- c(
  "((log(t) - ax) / bx)",
  "(exp(-(exp(zx))))",
  "((1 / (bx * t)) * exp(zx - exp(zx)))",
  "((1 / (bx * t)) * exp(zx))",
  "(exp(zx))"
)

# Private function, only for package
# Gumbel functions for; zx, Sx, fx, hx, Hx in terms of zx
P_Gumbel_functions_zx <- c(
  "((t - ax) / bx)",
  "(exp(-(exp(zx))))",
  "((1 / bx) * exp(zx - exp(zx)))",
  "((1 / bx) * exp(zx))",
  "(exp(zx))"
)

# Private function, only for package
# Frechet functions for; zx, Sx, fx, hx, Hx in terms of zx
P_Frechet_functions_zx <- c(
  "((log(t) - ax) / bx)",
  "(1 - exp(-(exp(-(zx)))))",
  "((1 / (bx * t)) * exp(-(zx) - exp(-(zx))))",
  "((1 / (bx * t)) * exp(-(zx) - exp(-(zx)))) / (1 - exp(-(exp(-(zx)))))",
  "-(log(1 - exp(-(exp(-(zx))))))"
)


# private functions converting functions above 'zx' to those in 'axbx'
P_Weibull_functions_axbx <- unlist(lapply(P_Weibull_functions_zx, function(x) gsub("zx", P_Weibull_functions_zx[1], x)))
P_Gumbel_functions_axbx <- unlist(lapply( P_Gumbel_functions_zx, function(x) gsub("zx",  P_Gumbel_functions_zx[1], x)))
P_Frechet_functions_axbx <- unlist(lapply(P_Frechet_functions_zx, function(x) gsub("zx", P_Frechet_functions_zx[1], x)))

# private function naming functions in list above
# name functions; 'z', 'S', 'f', 'h', 'H'
P_Weibull_z_S_f_h_H_axbx <- setNames(P_Weibull_functions_axbx, c("z", "S", "f", "h", "H"))
P_Gumbel_z_S_f_h_H_axbx <- setNames( P_Gumbel_functions_axbx, c("z", "S", "f", "h", "H"))
P_Frechet_z_S_f_h_H_axbx <- setNames(P_Frechet_functions_axbx, c("z", "S", "f", "h", "H"))

# private function making list of lists above
# create list of named functions
P_WGF_functions_axbx <- list(P_Weibull_z_S_f_h_H_axbx, P_Gumbel_z_S_f_h_H_axbx, P_Frechet_z_S_f_h_H_axbx)

# private function
# name list above according to probability distribution
P_WGF_named_function_list_z_S_f_h_H <- setNames(P_WGF_functions_axbx, c("Weibull", "Gumbel", "Frechet"))

# private functions
# call survival function in axbx by names of distribution & function
P_get_func_by_name_of_distribution_and_function <- function (dx = "", survfunc = ""){
  P_WGF_named_function_list_z_S_f_h_H[[dx]][survfunc]
}


# Private function, only for package
# writes S(t) for eval/paste/text into uniroot
# used in 1st version of etd_H1 model
P_eqn_Sx<- function(dx = ""){
  Sx <- "Sx"
  if((dx == "Weibull") | (dx == "weibull")){
    Sx <- "exp(-exp((log(t) - a1) / b1))"
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    Sx <- "exp(-exp((t - a1) / b1))"
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    Sx <- "1 - exp(-(exp(-(((log(t) - a1) / b1)))))"
  }  else {
    Sx <- "exp(-exp((log(t) - a1) / b1))"
  }
  return(Sx)
}


# Private function, only for package
# writes H(t) for eval/paste/text into uniroot
# used in 1st version of etd model
P_eqn_Hx <- function(dx = ""){
  Hx <- "Hx"
  if((dx == "Weibull") | (dx == "weibull")){
    Hx <- "-(log(exp(-exp((log(t) - a1) / b1))))"
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    Hx <- "-(log(exp(-exp((t - a1) / b1))))"
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    Hx <- "-(log(1 - exp(-(exp(-(((log(t) - a1) / b1)))))))"
  }  else {
    Hx <- "-(log(exp(-exp((log(t) - a1) / b1))))"
  }
  return(Hx)
}


# Private function, only for package
# Get zx as function of distribution type
# NB defaults to Weibull
P_get_zx_logscale <- function(t, log.ax, log.bx, dx){
  zx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    zx <- (log(t) - exp(log.ax)) / exp(log.bx)
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    zx <- (t - exp(log.ax)) / exp(log.bx)
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    zx <- (log(t) - exp(log.ax)) / exp(log.bx)
  } else {
    zx <- (log(t) - exp(log.ax)) / exp(log.bx)
  }
  return(zx)
}


# Private function, only for package
# Get fx as function of distribution type
P_get_fx_logscale <- function(t, zx, log.bx, dx){
  fx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    fx <- (1 / (exp(log.bx) * t)) * exp(zx - exp(zx))
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    fx <- (1 / exp(log.bx)) * exp(zx - exp(zx))
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    fx <- (1 / (exp(log.bx) * t)) * exp(-zx - exp(-zx))
  }  else {
    fx <- (1 / (exp(log.bx) * t)) * exp(zx - exp(zx))
  }
  return(fx)
}

# Private function, only for package
# Get hx as function of distribution type
P_get_hx_logscale <- function(t, zx, log.bx, dx){
  hx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    hx <- (1 / (exp(log.bx) * t)) * exp(zx)
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    hx <- (1 / exp(log.bx)) * exp(zx)
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    hx <- (1 / (exp(log.bx) * t)) * exp(-zx - exp(-zx)) / (1 - exp(-exp(-zx)))
  }  else {
    hx <- (1 / (exp(log.bx) * t)) * exp(zx)
  }
  return(hx)
}

# Private function, only for package
# Get Sx as function of distribution type
P_get_Sx_logscale <- function(t, zx, dx){
  Sx <- 0
  if((dx == "Weibull") | (dx == "weibull")){
    Sx <- exp(-exp(zx))
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    Sx <- exp(-exp(zx))
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") | (dx == "fr\u00E9chet")) {
    Sx <- 1 - exp(-exp(-zx))
  }  else {
    Sx <- exp(-exp(zx))
  }
  return(Sx)
}

### set defaults a1b1a2b2 ###
# Private function, only for package
# functions/expressions below used to set up default expressions
# for a1, b1, a2, b2

# empty functions used in nll function
P_func_a1b1a2b2 = list(
  func_a1 = function(){},
  func_b1 = function(){},
  func_a2 = function(){},
  func_b2 = function(){}
)

# default expressions
P_default_expressions_a1b1a2b2 = alist(
  exp_a1 = expression(a1),
  exp_b1 = expression(b1),
  exp_a2 = expression(a2),
  exp_b2 = expression(b2)
)

# expression for default function expressions
P_body_expressions_a1b1a2b2 <- expression(
  body(func_a1) <- eval(exp_a1),
  body(func_b1) <- eval(exp_b1),
  body(func_a2) <- eval(exp_a2),
  body(func_b2) <- eval(exp_b2)
)

# set formals as 'list of variables to estimate'
P_formal_expressions_a1b1a2b2 <- expression(
  formals(func_a1) <- formals(func_b1) <- formals(func_a2) <- formals(func_b2) <- default_list_of_variables_to_estimate_01
)


#############################

### set defaults a1b1a2b2p1 ###
# Private function, only for package
# functions/expressions below used to set up default expressions
# for a1, b1, a2, b2, p1

# empty functions used in nll function
P_func_a1b1a2b2p1 = list(
  func_a1 = function(){},
  func_b1 = function(){},
  func_a2 = function(){},
  func_b2 = function(){},
  func_p1 = function(){}
)

# default expressions
P_default_expressions_a1b1a2b2p1 = alist(
  exp_a1 = expression(a1),
  exp_b1 = expression(b1),
  exp_a2 = expression(a2),
  exp_b2 = expression(b2),
  exp_p1 = expression(p1)
)

# expression for default function expressions
P_body_expressions_a1b1a2b2p1 <- expression(
  body(func_a1) <- eval(exp_a1),
  body(func_b1) <- eval(exp_b1),
  body(func_a2) <- eval(exp_a2),
  body(func_b2) <- eval(exp_b2),
  body(func_p1) <- eval(exp_p1)
)

# set formals as 'list of variables to estimate'
P_formal_expressions_a1b1a2b2p1 <- expression(
  formals(func_a1) <- formals(func_b1) <- formals(func_a2) <- formals(func_b2) <- formals(func_p1) <- default_list_of_variables_to_estimate_02
)




#############################


#############################

### set defaults a1b1a2b2theta ###
# Private function, only for package
# functions/expressions below used to set up default expressions
# for a1, b1, a2, b2, theta

# empty functions used in nll function
P_func_a1b1a2b2theta = list(
  func_a1 = function(){},
  func_b1 = function(){},
  func_a2 = function(){},
  func_b2 = function(){},
  func_theta = function(){}
)

# default expressions
P_default_expressions_a1b1a2b2theta = alist(
  exp_a1 = expression(a1),
  exp_b1 = expression(b1),
  exp_a2 = expression(a2),
  exp_b2 = expression(b2),
  exp_theta = expression(theta)
)

# expression for default function expressions
P_body_expressions_a1b1a2b2theta <- expression(
  body(func_a1) <- eval(exp_a1),
  body(func_b1) <- eval(exp_b1),
  body(func_a2) <- eval(exp_a2),
  body(func_b2) <- eval(exp_b2),
  body(func_theta) <- eval(exp_theta)
)

# set formals as 'list of variables to estimate'
P_formal_expressions_a1b1a2b2theta <- expression(
  formals(func_a1) <- formals(func_b1) <- formals(func_a2) <- formals(func_b2) <- formals(func_theta) <- default_list_of_variables_to_estimate_03
)


#############################

# private function used in estimating confidence intervals
# used in function estimating confidence intervals
P_get_dh_da <- function(t, a2, b2, z2, d2){

  dh_da <- 0

  if((d2 == "Weibull") | (d2 == "weibull")){

    dh_da <- -(1 / ((b2^2) * t)) * exp(z2)

  } else if ((d2 == "Gumbel") | (d2 == "gumbel")) {

    dh_da <- -(1 / b2^2) * exp(z2)

  } else if ((d2 == "Frechet") | (d2 == "frechet") | (d2 == "Fr\u00E9chet") | (d2 == "fr\u00E9chet")) {

    A <- exp(-z2 - exp(-z2)) * ((1 / b2) - (exp(-z2) / b2))
    B <- (b2 * t) * (1 - exp(-exp(-z2)))
    C <- exp((-2 * z2) - 2 * exp(-z2))
    D <- ((b2^2 * t) * ((1 - exp(-exp(-z2)))^2))

    dh_da <- (A / B) - (C / D)

    #   dh_da <- (((1 / b2) - (exp(-z2) / b2)) * exp(-z2 - exp(-z2))) / ((b2 * t)*(1 - exp(-exp(-z2)))) - ((exp(-2 * z2 - 2 * exp(-z2))) / (b2^2 * t * (1 - exp(-exp(-z2)))))

  } else {
    stop ("No probability distribution given for d2")
  }
  return(dh_da)
}

# private function used in estimating confidence intervals
# used in function estimating confidence intervals
P_get_dh_db <- function(t, a2, b2, z2, d2){

  dh_db <- 0

  if((d2 == "Weibull") | (d2 == "weibull")){

    dh_db <- (-1 / ((b2^2) * t)) * (z2 * exp(z2) - exp(z2))

  } else if ((d2 == "Gumbel") | (d2 == "gumbel")) {

    dh_db <- (-1 / b2^2) * exp(z2) - ((t - a2) * exp(z2) / b2^3)

  } else if ((d2 == "Frechet") | (d2 == "frechet") | (d2 == "Fr\u00E9chet") | (d2 == "fr\u00E9chet")) {

    A <- (exp(-z2 - exp(-z2)) * ((log(t) - a2) / b2^2) - (((log(t) - a2) * exp(-z2)) / b2^2))
    B <- ((b2 * t) * (1 - exp(-exp(-z2))))
    C <- exp(-z2 - exp(-z2))
    D <- ((b2^2 * t) * (1 - exp(-exp(-z2))))
    E <- (log(t) - a2) * exp(-2 * z2 - 2 * exp(-z2))
    G <- ((b2^3 * t) * (1 - exp(-exp(-z2)))^2)

    dh_db <- (A / B) - (C / D) - (E / G)

  } else {
    stop ("No probability distribution given for d2")
  }
  return(dh_db)
}


# private function
# used in generating random times of death
P_generator_random_times <- function(ax, bx, dx = ""){

  x <- 0
  r <- stats::runif(1)

  if((dx == "Weibull") | (dx == "weibull")){
    x <- exp((bx * log(-log(r))) + ax)
  } else if ((dx == "Gumbel") | (dx == "gumbel")) {
    x <- bx * log(-log(r)) + ax
  } else if ((dx == "Frechet") | (dx == "frechet") | (dx == "Fr\u00E9chet") |
             (dx == "fr\u00E9chet")) {
    x <- exp((bx * log(-log(1 - r))) + ax)
  } else {
    stop('probability distribution not defined')
  }
  return(x)
}


# Private function used by 'sim_data_nll_basic'
# function sequentially checks input values are numeric and greater than zero
# stops & returns warning when not
P_check_inputs_numeric_and_greater_than_zero <- function(...){

  input <- c(as.list(environment()), list(...))
  input_length <- length(input)
  input_names <- names(input)

  i <- 1

  for(i in 1:input_length){

    if(is.numeric(input[[i]]) == FALSE) stop(paste(input_names[i], 'must be numeric'), call. = FALSE)
    if(input[[i]] <= 0) stop(paste(input_names[i], 'must be > 0'), call. = FALSE)
  }
}


# Private function, only for package
# Gets expression for hazard function
# in terms of a2, b2, d2
# used by 'anovir::conf_ints_virulence'

P_hazard_expression <- function(d2 = ""){
  if((d2 == "Weibull") | (d2 == "weibull")){
    expression(1 / (b2*t) * exp((log(t) - a2) / b2))
  } else if((d2 == "Gumbel") | (d2 == "gumbel")){
    expression((1 / b2) * exp((t - a2) / b2))
  } else if((d2 == "Frechet") | (d2 == "frechet") |
            (d2 == "Fr\u00E9chet") | (d2 == "fr\u00E9chet")) {
    expression((1 / (b2 * t)) * exp(-((log(t) - a2) / b2) - exp(-((log(t) - a2) / b2))) / (1 - exp(-exp(-((log(t) - a2) / b2)))))
  } else {
    stop("Please specify probability distribution d2: Weibull, Gumbel or Frechet")
  }
}

# P_hazard_expression(d2 = 'Weibull')
# P_hazard_expression(d2 = 'Gumbel')
# P_hazard_expression(d2 = 'Frechet')
# P_hazard_expression('Weibull')
# P_hazard_expression()
# P_hazard_expression(d2 = 'mango')



