kl_helpers_surv <- function(fam) {
  # NOTE: we should get rid off these, they are not much of a help..
  kl_surv <- function(pref, data, psub)
    colSums( abs( (exp(psub$mu - pref$mu) - 1) / (exp(psub$mu - pref$mu) + 1) ) )

  dis_na <- function(pref, psub, wobs) rep(0, ncol(pref$mu)) 
  predvar_na <- function(mu, dis, wsample=1) { 0 }
  
  ll_surv <- function(mu, dis, y, weights=1, ll_args_surv) {
    time <- y[,1]
    status <- y[,2]
    
    args <- list(basehaz   = ll_args_surv$basehaz,
                 aux       = ll_args_surv$aux,
                 intercept = ll_args_surv$alpha,
                 time      = time )
    mu <- t(mu) - c(ll_args_surv$alpha)
    lhaz  <- do.call(rstanarm:::evaluate_log_basehaz,  args) + mu
    lsurv <- do.call(rstanarm:::evaluate_log_basesurv, args) * exp(mu)
    l_lik <- lsurv
    l_lik[ ,status == 1] <- lhaz[ ,status == 1] + lsurv[ ,status == 1] 
    t(as.matrix(weights*l_lik))
  }
  
  dev_surv <- function(mu, dis, y, weights=1, ll_args_surv) {
    time <- y[,1]
    status <- y[,2]
    
    args <- list(basehaz   = ll_args_surv$basehaz,
                 aux       = ll_args_surv$aux,
                 intercept = ll_args_surv$alpha,
                 time      = time )
    mu <- t(mu) - c(ll_args_surv$alpha)
    lhaz  <- do.call(rstanarm:::evaluate_log_basehaz,  args) + mu
    lsurv <- do.call(rstanarm:::evaluate_log_basesurv, args) * exp(mu)
    l_lik <- lsurv
    l_lik[ ,status == 1] <- lhaz[ ,status == 1] + lsurv[ ,status == 1]    -2*weights*l_lik
  }
  
  ppd_weibull <- function(mu, dis, weights = 1, alpha, aux, basehaz ) {
    rweibull(length(mu), aux, exp( - mu / aux ))
    # sim_tte(mu, alpha, aux,  basehaz, obs.only = F,
    #         obs.aug = T,
    #         delta = 0.01,
    #         time = 200,
    #         end_time = 200)
  }
  
  # function for computing mu = E(y)
  mu_fun <- function(x, alpha, beta, offset) {
    if (!is.matrix(x)) stop('x must be a matrix.')
    if (!is.matrix(beta)) stop('beta must be a matrix')
    fam$linkinv(cbind(1, x) %*% rbind(alpha, beta) + offset)
  }
  
  latent_fun <- function(x, alpha, beta, offset) {
    if (!is.matrix(x)) stop('x must be a matrix.')
    if (!is.matrix(beta)) stop('beta must be a matrix')
    return(cbind(1, x) %*% rbind(alpha, beta) + offset)
  }
  c(list(kl = kl_surv, ll_fun = ll_surv, deviance = dev_surv, dis_fun = dis_na, predvar = predvar_na, ppd_fun = ppd_weibull),
  list(mu_fun = mu_fun, latent_fun = latent_fun), fam)
}

surv_family <- function(){
  out <- list(
    family = "surv",
    link   = "extreme",
    linkfun = function(mu) log(mu),
    linkinv = function(eta) exp(eta),
    variance = function() NULL,
    ## The following feels are not completed but are not used in the package
    dev.resids = function() NULL,
    aic      = function() NULL,
    mu.eta   = function(eta) NULL, 
    validmu  = function(mu) TRUE,
    valideta  = function(mu) TRUE,
    simulate  = NULL
  )
  structure(out, class = "family")
}


is_stan_surv <- function(fit) {
  "stansurv" %in% class(fit)
}

is_surv_family <- function(x){
  if(is.character(x$family)){
    grepl("surv", x$family)
  } else {
    grepl("surv", x$family$family)
  }
  
}

posterior_survlinpred <- function(object, newdata = NULL, transform=T, offset ){
  if(is.null(newdata)) {
    zt <- object$x
  } else {
    zt <- newdata
  }
  draws <- rstanarm:::ll_args.stansurv(object)$draws
  alpha <- draws$alpha
  beta <- draws$beta
  res <-  cbind(1, zt) %*% t(cbind(alpha, beta))
  exp(res)
}



posterior_survtimepred <- function(object, newdata = NULL, transform=T, offset ){
  if(is.null(newdata)) {
    zt <- object$x
  } else {
    zt <- newdata
  }
  draws <- rstanarm:::ll_args.stansurv(object)$draws
  alpha <- draws$alpha
  beta <- draws$beta
  shape <- draws$aux
  res <-exp( - ( cbind(1, zt) %*% t(cbind(alpha, beta))  ) / shape  )
  res
}


pseudo_data_surv <- function(mu, family, shape, weights = wobs) {
  if(is.null(shape)) {
    ## exponential model
    shape = 1
  }
  t( exp( - t(mu) / shape ) )
}



# project_survival <-  function(vind, p_ref, d_train, family_kl, intercept, regul=1e-9, coef_init=NULL) {
#   
#   mu <- p_ref$mu_latent
#   dis <- p_ref$dis
#   aux <- exp(p_ref$log_aux)
#   
#   if (is.null(refmodel$wobs)) {
#     wobs <- rep(1.0, NROW(mu))
#   } else {
#     wobs <- refmodel$wobs
#   }
#   
#   if (is.null(p_ref$weights)) {
#     wsample <- rep(1.0, NCOL(mu))
#   } else {
#     wsample <- p_ref$weights
#   }
#   
#   wobs <- wobs / sum(wobs)
#   wsample <- wsample / sum(wsample)
#   
#   form <- as.formula(paste0( "Surv(time, status) ~ ", paste( paste0("V", seq_along(vind), collapse = " + " ) ) ) ) 
#   
#   pobs <- pseudo_data_surv(mu, family, aux, weights = wobs)
#   x_train <- d_train$z[, vind, drop = FALSE]
#   status <- matrix(1, nrow(pobs), ncol(pobs))
#   
#   ## should add censoring, or is losing information 
#   # max_time_train <- max(d_train$y[,1])
#   # status[pobs > max_time_train] <- 0
#   # pobs[pobs > max_time_train] <- max_time_train
#   
#   proj_refit <- lapply(seq_len(ncol(pobs)), function(j) {
#     survdata <- cbind.data.frame( data.frame( time = pobs[,j], status = status[,j]), x_train)
#     survival::survreg(form, survdata)
#   } )
#   
#   shape <- sapply(proj_refit, function(x) 1/x$scale )
#   coefs <- sapply(proj_refit, function(x) - x$coefficients/x$scale)
#   
# }

