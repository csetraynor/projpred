#' @export
surv_family <- function(fit){
  basehaz_ok <- c("exponential", "weibull")
  basehaz_name <- rstanarm:::get_basehaz_name(fit)
  if(basehaz_name %notin% basehaz_ok) 
    stop("Baehaz must be exponential or weibull.")
  surv_family <- switch(
    basehaz_name,
    "exponential" = surv_exp_family(),
    "weibull"     = surv_wei_family()
    )
  return(surv_family)
}

surv_exp_family <- function(){
  out <- list(
    family = "surv_exponential",
    link   = "extreme",
    linkfun = function(mu) exp(mu),
    linkinv = function(eta) log(eta),
    variance = function() NULL,
    dev.resids = function() NULL,
    aic      = function() NULL,
    mu.eta   = function(eta) rep.int(1, length(eta)),
    validmu  = function(mu) TRUE,
    valideta  = function(mu) TRUE,
    simulate  = NULL
  )
  structure(out, class = "family")
}

surv_wei_family <- function(){
  out <- list(
    family = "surv_weibull",
    link   = "extreme",
    linkfun = function(mu) mu,
    linkinv = function(mu) mu,
    variance = function() NULL,
    dev.resids = function() NULL,
    aic      = function() NULL,
    mu.eta   = function(eta) rep.int(1, length(eta)),
    validmu  = function(mu) TRUE,
    valideta  = function(mu) TRUE,
    simulate  = NULL
  )
  structure(out, class = "family")
}

is_stan_surv <- function(fit) {
  "stansurv" %in% class(fit)
}

post_surv_latentfactor <- function(fit, newdata){
  args <- rstanarm:::ll_args.stansurv(fit, newdata = newdata)
  x <- rstanarm:::.xdata_surv(args$data) 
  x <- as.matrix(x)
  beta <- args$draws$beta
  x %*% t(beta)
}

post_surv_exponential_timepred <- function(fit, newdata) {
  args <- rstanarm:::ll_args.stansurv(fit, newdata = newdata)
  x <- rstanarm:::.xdata_surv(args$data) 
  x <- cbind(rep(1, nrow(x)), as.matrix(x))
  beta <- cbind(args$draws$alpha, args$draws$beta)
  mu <-   x %*% t(beta)
  exp ( - mu )
}

get_shape <- function(fit){
  args <- rstanarm:::ll_args.stansurv(fit)
  shape <- args$draws$aux
}

post_surv_exponential_timepred <- function(fit, newdata) {
  args <- rstanarm:::ll_args.stansurv(fit, newdata = newdata)
  x <- rstanarm:::.xdata_surv(args$data) 
  x <- cbind(rep(1, nrow(x)), as.matrix(x))
  beta <- cbind(args$draws$alpha, args$draws$beta)
  mu <-   x %*% t(beta)
  shape <- c(args$draws$aux)
  exp ( - mu %*% diag(1 / shape) )
}