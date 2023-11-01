require(dplyr)
require(TMB)

# Read in data function
read.data <- function(filename) {
  dat <- readLines(filename)
  ls.out <- list()
  
  ls.out$ncreeks <- gsub('\t', '', dat[2]) %>% as.numeric()
  ls.out$nrel <- gsub('\t', '', dat[4])  %>% as.numeric()
  ls.out$nsites <- gsub('\t', '', dat[6]) %>% as.numeric()
  ls.out$nperiods <- gsub('\t', '', dat[8]) %>% as.numeric()
  ls.out$ntraps <- gsub('\t', '', dat[10]) %>% as.numeric()
  
  count.ls <- with(ls.out,
                   strsplit(dat[11 + 1:(ncreeks*nsites*nperiods)], '\t'))
  ls.out$count.mat <- do.call(rbind, count.ls)[,1:ls.out$ntraps] %>% apply(2, as.numeric)
  
  ls.out$distances <- with(ls.out,
                           strsplit(dat[13 + ncreeks*nsites*nperiods], '\t')[[1]] %>% 
                             as.numeric())
  ls.out$times <- with(ls.out,
                       strsplit(dat[15 + ncreeks*nsites*nperiods], '\t')[[1]] %>% 
                         as.numeric())
  return(ls.out)
}

# Set up model function
# dat.ls <- dat; disp.model <- disp.mod; count.model <- 'neg.binom'; dist.cutoff <- 50; sigma.type <- "constant"
make.inputs <- function(dat.ls, disp.model, count.model, dist.cutoff, sigma.type) {
  disp.model.num <- switch(disp.model,
                           normal = 1,
                           exponential = 2,
                           cauchy = 3,
                           no.dispersal = 4)
  if(is.null(disp.model.num)) 
    stop('Please choose "normal", "exponential", "cauchy", or "no.dispersal" for dispersal model')
  
  count.model.num <- switch(count.model,
                            poisson = 1,
                            neg.binom = 2)
  if(is.null(count.model.num))
    stop('Please choose "poisson" or "neg.binom" for count model')
  
  if(!(sigma.type %in% c('random', 'fixed', 'asympBH', 'asympVB', 'constant', 'none')))
    stop('Please choose "random", "fixed", "asympBH", "asympVB", "none", or "constant" for sigma type')
  
  Data <- with(dat.ls, 
               list(disp_model = disp.model.num,
                    count_model = count.model.num,
                    is_random = as.numeric(sigma.type == 'random'),
                    is_asymptotic = switch(sigma.type, asympBH=1, asympVB=2, 0),
                    nrel = nrel,
                    nsites = nsites,
                    nperiods = nperiods,
                    ntraps = ntraps,
                    countmat = count.mat,
                    distances = distances,
                    times = times,
                    dist_cutoff = dist.cutoff))
  
  Map <- list()
  Random <- NULL
  
  if(count.model.num == 1) {
    Map$log_overdispersion <- factor(NA)
  }
  
  if(disp.model=='no.dispersal') {
    Map$log_sig_disp_mu <- factor(NA) 
    Map$log_sig_disp_sig <- factor(NA)
    Map$log_sig_disp_eps <- rep(factor(NA), dat.ls$nperiods)
  }
  if(sigma.type == 'random') {
    Random <- 'log_sig_disp_eps'
  } else {
    if(sigma.type == 'fixed') {
      Map$log_sig_disp_mu <- factor(NA) 
      Map$log_sig_disp_sig <- factor(NA)
    }
    if(sigma.type == 'constant') {
      Map$log_sig_disp_eps <- rep(factor(NA), dat.ls$nperiods)
      Map$log_sig_disp_sig <- factor(NA)
    }
    if(sigma.type == 'asympBH' | sigma.type == 'asympVB') {
      Map$log_sig_disp_eps <- rep(factor(NA), dat.ls$nperiods) 
    }
  }
  
  
  return(list(Data = Data, Map = Map, Random = Random))
}

# Compile and load model
#compile('DM_const_sig.cpp')
try(dyn.unload(dynlib("new_model/code/helper_functions/DM_MM_sig")))
compile('new_model/code/helper_functions/DM_MM_sig.cpp')

#dyn.load(dynlib("DM_const_sig"))
dyn.load(dynlib("new_model/code/helper_functions/DM_MM_sig"))

Parameters <- list(
  survival = .9,
  detectability = .7,
  sig_disp = 20,
  overdispersion = 1
)


