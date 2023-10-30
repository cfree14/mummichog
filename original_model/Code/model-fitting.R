# Fit models to data ------------------------------------------------------

source("original_model/Code/for_TMB.R")

aic.emp <- array(0, dim=c(4,3,4), dimnames = list(creek=NULL, 
                                              disp.mod=c('normal', 'exponential', 'cauchy'),
                                              mod = c('constant', 'fixed effect', 'random effect', 'asymptote')))
emp.mod.fits <- list()

aic.no.dispersal <- numeric(4)
no.dispersal.fits <- list()

# Fit models
for(creek in 1:4) {
  filename <- paste('original_model/data/crk', creek, '.dat', sep='')
  dat <- read.data(filename)
  
  Parameters <- list(
    logit_survival = log(.98/.02),
    log_detectability = log(1),
    log_sig_disp_mu = log(20),
    log_sig_disp_sig = log(1),
    log_sig_disp_eps = rep(log(1), dat$nperiods), #log(rep(20, 9)),
    log_overdispersion = log(2)
  )
  
  emp.mod.fits[[creek]] <- list()
  
  for(disp.mod in c('normal', 'exponential', 'cauchy')) {
    emp.mod.fits[[creek]][[disp.mod]] <- list()
    
    for(mod in 1:4) {
      to.fit <- make.inputs(dat, disp.mod, 'neg.binom',
                            50, switch(mod, 'constant', 'fixed', 'random', 'asympBH'))
      Parameters$log_sig_disp_mu <- log(30)
      Parameters$log_sig_disp_eps <- rep(log(30), dat$nperiods)
      if(mod == 1) Parameters$log_sig_disp_eps <- rep(0, dat$nperiods)
      if(mod == 2) Parameters$log_sig_disp_mu <- 0
      
      model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_MM_sig", map=to.fit$Map,
                         random = to.fit$Random)
      model$env$beSilent()
      Opt <- tryCatch(nlminb(start=model$par, objective=model$fn, gradient=model$gr),
                      error = function(e) NA)
      aic.emp[creek, disp.mod, mod] <- ifelse(!is.na(Opt[1]), 
                                              ifelse(sdreport(model)$pdHess, TMBhelper::TMBAIC(Opt),
                                              NA), NA)
      emp.mod.fits[[creek]][[disp.mod]][[mod]] <- model
    }
  }
  to.fit <- make.inputs(dat, 'no.dispersal', 'neg.binom',
                        50, 'none')
  Parameters$log_sig_disp_mu <- 0
  Parameters$log_sig_disp_eps <- rep(0, dat$nperiods)
  
  model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_MM_sig", map=to.fit$Map,
                     random = to.fit$Random)
  model$env$beSilent()
  Opt <- tryCatch(nlminb(start=model$par, objective=model$fn, gradient=model$gr),
                  error = function(e) NA)
  aic.no.dispersal[creek] <- ifelse(!is.na(Opt[1]), 
                                    ifelse(sdreport(model)$pdHess, TMBhelper::TMBAIC(Opt),
                                           NA), NA)
  no.dispersal.fits[[creek]] <- model
}
