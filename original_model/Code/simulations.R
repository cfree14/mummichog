require(TMB)
require(MASS)
set.seed(789890235)
dyn.load(dynlib("Code/DM_MM_sig"))

# Data structure
dat <- read.data('Data/crk2.dat')
ncreeks <- dat$ncreeks
nrel <- dat$nrel
nsites <- dat$nsites
nperiods <- dat$nperiods
ntraps <- dat$ntraps
distances <- dat$distances
times <- dat$times

# True values
survival <- 0.99
detectability <- 2
sig.disp <- 32.5*times/(times+6.35)
#rexp(nperiods, .2) %>% sort(decreasing = TRUE) %>% cumsum() %>% sort() #rep(20, nperiods) #
annual.mort <- (-log(survival)*365)
t.max <-  exp((1.46-log(annual.mort))/1.01)
# t.max <- 4.22/(-log(survival)*365)
overdispersion <- 2

disp.mods <- c('normal', 'exponential', 'cauchy')
disp.structures <- c('constant', 'fixed', 'random', 'asympBH', 'asympVB')

nreps <- 500
nmods <- 3
count.mat.sim <- array(0, dim = c(nmods, nsites*nperiods, ntraps),
                       dimnames = list(mod=NULL, site.period = NULL, trap=NULL))

aic.sim <- array(0, dim=c(nmods, nmods, nreps), 
                 dimnames=list(sim.mod=disp.mods, est.mod=disp.mods, 
                               rep=1:nreps))

fitted.mods <- list()
for(sim.mod in 1:nmods) {
  fitted.mods[[sim.mod]] <- list()
  for(est.mod in 1:(nmods+1)) {
    fitted.mods[[sim.mod]][[est.mod]] <- list()
    if(est.mod != 4) {
      for(struc in 1:length(disp.structures)) {
        fitted.mods[[sim.mod]][[est.mod]][[struc]] <- list()
      }
    }
  }
}

simulate.counts <- function(distance, times, disp.mod, ntraps, nrel, survival, detectability,
                            overdispersion = NA, half.distn, ...) {
  ## Take dispersal model choice, other parameters, to simulate observed counts at each of
  ## ntraps traps. Returns vector of length ntraps.
  ## If overdispersion parameter is numeric, simulates negative binomials counts. 
  ## Otherwise (NA, NULL, FALSE, etc.) simulates Poison counts.
  dist.factor <- (1+half.distn) * disp.mod(distance, ...)
  pred.count <- nrel * survival^times * detectability * dist.factor
  if(is.numeric(overdispersion)) {
    counts <- rnbinom(ntraps, mu = pred.count, size = overdispersion)
  } else {
    counts <- rpois(ntraps, pred.count)
  }
  return(counts)
}

Parameters <- list(
  logit_survival = log(.99/.01),
  log_detectability = log(2),
  log_sig_disp_mu = log(30),
  log_sig_disp_sig = log(1),
  log_sig_disp_eps = rep(0, 9),
  log_overdispersion = log(2)
)



for(ii in 1:nreps) {
  # simulate data
  ind <- 1
  for(period in 1:nperiods) {
    for(site in 1:nsites) {
      # half-normal
      count.mat.sim[1,ind,] <- simulate.counts(distance = distances[site], 
                                               times=times[period], disp.mod = dnorm,
                                               ntraps = ntraps, nrel = nrel,
                                               survival = survival, 
                                               detectability = detectability, 
                                               overdispersion = overdispersion,
                                               half.distn = TRUE, mean = 0, 
                                               sd = sig.disp[period])
      # exponential
      count.mat.sim[2,ind,] <- simulate.counts(distance = distances[site], 
                                               times=times[period], disp.mod = dexp,
                                               ntraps = ntraps, nrel = nrel,
                                               survival = survival, 
                                               detectability = detectability, 
                                               overdispersion = overdispersion,
                                               half.distn = FALSE,
                                               rate = 1/sig.disp[period])
      # half-cauchy
      count.mat.sim[3,ind,] <- simulate.counts(distance = distances[site], 
                                               times=times[period], disp.mod = dcauchy,
                                               ntraps = ntraps, nrel = nrel,
                                               survival = survival, 
                                               detectability = detectability, 
                                               overdispersion = overdispersion,
                                               half.distn = TRUE, location = 0, 
                                               scale = sig.disp[period])
      
      ind <- ind + 1
    }
  }
  # fit data
  for(sim.mod in 1:nmods) {
    for(est.mod in 1:nmods) {
      for(struc in 1:length(disp.structures)) {
        input.ls <- make.inputs(dat.ls = list(nrel = nrel,
                                              nsites = nsites,
                                              nperiods = nperiods,
                                              ntraps = ntraps,
                                              count.mat = count.mat.sim[sim.mod,,],
                                              distances = distances,
                                              times = times), 
                                disp.model = disp.mods[est.mod],
                                count.model = 'neg.binom', dist.cutoff = 50, 
                                sigma.type = disp.structures[struc]) 
        Parameters$log_sig_disp_mu <- log(20)
        Parameters$log_sig_disp_eps <- rep(log(20), nperiods)
        if(disp.structures[struc] == 'constant') Parameters$log_sig_disp_eps <- rep(0, nperiods)
        if(disp.structures[struc] == 'fixed') Parameters$log_sig_disp_mu <- 0
        
        model <- MakeADFun(data = input.ls$Data, parameters = Parameters, 
                           map = input.ls$Map, 
                           DLL="DM_MM_sig")
        model$env$beSilent()
        Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr)
        fitted.mods[[sim.mod]][[est.mod]][[struc]][[ii]] <- model
        if(struc == (which(disp.structures=='asympBH')))
          aic.sim[sim.mod, est.mod, ii] <- TMBhelper::TMBAIC(Opt)
      }
    }
      recaps <- apply(count.mat.sim[sim.mod,,], 1, sum) %>%
        tapply(rep(1:nperiods, each=nsites), sum) %>% as.vector()
      
      # fit dispersal free model
      input.ls <- make.inputs(dat.ls = list(nrel = nrel,
                                            nsites = nsites,
                                            nperiods = nperiods,
                                            ntraps = ntraps,
                                            count.mat = count.mat.sim[sim.mod,,],
                                            distances = distances,
                                            times = times), 
                              disp.model = 'no.dispersal',
                              count.model = 'neg.binom', dist.cutoff = 50, 
                              sigma.type = 'none') 
      Parameters$log_sig_disp_mu <- 0
      Parameters$log_sig_disp_eps <- rep(0, nperiods)
      
      model <- MakeADFun(data = input.ls$Data, parameters = Parameters, 
                         map = input.ls$Map, 
                         DLL="DM_MM_sig")
      model$env$beSilent()
      Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr)
      fitted.mods[[sim.mod]][[nmods+1]][[ii]] <- model
    
  }
}

true.val.mat <- matrix(c(qnorm(.75, 0, sig.disp[1]), 
                         2*(pnorm(50, 0, sig.disp[1]) - 0.5),
                         qexp(.5, 1/sig.disp[1]), 
                         pexp(50, 1/sig.disp[1]),
                         qcauchy(.75, 0, sig.disp[1]), 
                         2*(pcauchy(50, 0, sig.disp[1]) - 0.5)),
                       nrow=3, byrow=TRUE,
                       dimnames = list(model = disp.mods,
                                       val = c('fifty_pct', 'pct_at_dist')))

# Collect survival and dispersal paramters. All models. 
sim.res <- sim.ls <- list()
surv.mat <- fifty.pct.mat <- pct.at.dist.mat <- matrix(0, nrow=nreps, ncol=length(disp.structures),
                   dimnames = list(rep=NULL, disp.str=disp.structures))

for(sim.mod in 1:nmods) {
  sim.res[[sim.mod]] <- list()
  for(est.mod in 1:nmods) {
    for(struc in 1:length(disp.structures)) {
      sdreports <- sapply(fitted.mods[[sim.mod]][[est.mod]][[struc]], sdreport)
      surv.mat[,struc] <- apply(sdreports, 2, function(x) 
        if(x$pdHess) {
          x['value']$value['annual_mort']
        } else NA)
      pct.at.dist.mat[,struc] <- apply(sdreports, 2, function(x) 
        if(x$pdHess) {
          x['value']$value['pct_at_dist']
        } else NA)
      fifty.pct.mat[,struc] <- apply(sdreports, 2, function(x) 
        if(x$pdHess) {
          x['value']$value['fifty_pct']
        } else NA)
    }
    sim.res[[sim.mod]][[est.mod]] <- data.frame(annual.mort=as.vector(surv.mat),
                                                 pct.at.dist=as.vector(pct.at.dist.mat),
                                                 fifty.pct=as.vector(fifty.pct.mat),
                                                 disp.str=rep(disp.structures, each=nreps),
                                                 est.mod=disp.mods[est.mod],
                                                 sim.mod=disp.mods[sim.mod],
                                                 stringsAsFactors = FALSE)
  }
  
  sdreports <- sapply(fitted.mods[[sim.mod]][[nmods+1]], sdreport)
  surv.vec <- apply(sdreports, 2, function(x) 
    if(x$pdHess) {
      x['value']$value['annual_mort']
    } else NA)
  sim.res[[sim.mod]][[nmods+1]] <- data.frame(annual.mort=surv.vec,
                                              pct.at.dist=NA,
                                              fifty.pct=NA,
                                              disp.str=NA,
                                              est.mod='no.dispersal',
                                              sim.mod=disp.mods[sim.mod],
                                              stringsAsFactors = FALSE)

  # from when no dispersal model was fit with glm.nb
  # sapply(fitted.mods[[sim.mod]][[nmods+1]],
  #                               function(x) 
  #                                 # tryCatch(-4.22/(coef(x)['times']*365),
  #                                 # tryCatch(exp((1.46-log(-(coef(x)['times']*365)))/1.01),
  #                                 tryCatch(-coef(x)['times']*365,
  #                                          error=function(e) return(NA))) %>%
  #   data.frame(annual.mort=., pct.at.dist=NA, fifty.pct=NA, disp.str=NA, est.mod='no dispersal', 
  #              sim.mod=disp.mods[sim.mod], stringsAsFactors = FALSE)
  #                                 # S = exp(nb glm estimate)
  #                                 # M = -log(S)
  #                                 # convert daily to annual rate
  #                                 # convert to max age
  #                                 # matrix(ncol=3, byrow=TRUE)
  sim.ls[[sim.mod]] <- do.call(rbind, sim.res[[sim.mod]][1:nmods]) %>%
    rbind(sim.res[[sim.mod]][[nmods+1]])
      # aic.sim[sim.mod, 4,] <- sapply(fitted.mods[[sim.mod]][[4]], function(mod)
  #   if(length(mod)>1) return(AIC(mod)) else return(NA))
  # 
  # temp[[sim.mod]] <- do.call(rbind, res[[sim.mod]]) %>% data.frame() %>% 
  #   mutate(est.mod = rep(c(disp.mods, 'no dispersal'), each=nreps))
}
sim.res.df <- do.call(rbind, sim.ls)
sim.res.df$disp.str[is.na(sim.res.df$disp.str)] <- 'no dispersal'

# save results ------------------------------------------------------------

save(sim.res.df, aic.sim, annual.mort, sig.disp, nreps, file = 'Code/sim-results.RData')
