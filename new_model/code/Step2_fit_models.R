

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(TMB)
library(TMBhelper)
library(tidyverse)
library(ggsidekick)
library(lubridate)

# Directories
codedir <- "new_model/code/helper_functions"
datadir <- "new_model/data"
plotdir <- "new_model/figures"

# Source helper functions
source(file.path(codedir, "for_TMB.R"))

# Fit to NJ data?
# FALSE - fits to Gulf of Mexico data
# TRUE - fits to NJ data
nj_yn <- T
if(nj_yn){
  # Creeks
  creeks <- c("Raritan-Downstream", "Raritan-Upstream", "Passaic-Downstream")
  files <- paste0(tolower(creeks), ".dat")
}else{
  # Creeks
  creeks <- paste("Creek", 1:4)
  files <- paste0("crk", 1:4, ".dat")
}


# Fit models
################################################################################

# Number of creeks
ncreeks <- length(creeks)

# Setup empty AIC container
# 3 creeks, 3 dispersal models, 4 model structures
aic.emp <- array(0, dim=c(ncreeks,3,4), dimnames = list(creek=NULL, 
                                              disp.mod=c('normal', 'exponential', 'cauchy'),
                                              mod = c('constant', 'fixed effect', 'random effect', 'asymptote')))
# Setup empty model fits constainer
emp.mod.fits <- list()

# Setup containers for "no dispersal" models 
aic.no.dispersal <- numeric(ncreeks)
no.dispersal.fits <- list()

# Loop through creeks and fit models
i <- 1
for(i in 1:ncreeks) {
  
  # Data structure
  # Datafile is a list with following:
  # ncreeks = # of creeks
  # nrel = # of released fish
  # nsites = # of sites
  # nperiods = # of time periods
  # ntraps = # of traps
  # count.mat = # of tagged fish at a date, site, trap
  # distances = distances of sites from tagging location (vector length of nsites)
  # times = days since time of tagging (vector length of nperiods)
  
  # Creek
  creek_do <- creeks[i]
  file_do <- files[i]
  
  # Read data
  dat <- read.data(file.path(datadir, file_do))
  
  # Setup initial parameters
  Parameters <- list(
    logit_survival = log(.98/.02),
    log_detectability = log(1),
    log_sig_disp_mu = log(20),
    log_sig_disp_sig = log(1),
    log_sig_disp_eps = rep(log(1), dat$nperiods), #log(rep(20, 9)),
    log_overdispersion = log(2)
  )
  
  # Setup empty model fit
  emp.mod.fits[[i]] <- list()
  
  # Loop through dispersal types
  # Half normal, exponential, half cauchy
  disp.mod <- "normal"
  for(disp.mod in c('normal', 'exponential', 'cauchy')) {
    
    # Setup empty model fit
    emp.mod.fits[[i]][[disp.mod]] <- list()
    
    # Loop through model structures
    # independent, random effect, 2-parameter asymptotic function
    mod <- 1
    for(mod in 1:4) {
      
      # Build data to fit model to
      to.fit <- make.inputs(dat, disp.mod, 'neg.binom',
                            50, switch(mod, 'constant', 'fixed', 'random', 'asympBH'))
      
      # Modify initial parameters
      Parameters$log_sig_disp_mu <- log(30)
      Parameters$log_sig_disp_eps <- rep(log(30), dat$nperiods)
      if(mod == 1) Parameters$log_sig_disp_eps <- rep(0, dat$nperiods)
      if(mod == 2) Parameters$log_sig_disp_mu <- 0
      
      # Fit model
      model <- TMB::MakeADFun(data=to.fit$Data, 
                              parameters=Parameters, 
                              DLL="DM_MM_sig", 
                              map=to.fit$Map,
                              random = to.fit$Random)
      model$env$beSilent()
      Opt <- tryCatch(nlminb(start=model$par, objective=model$fn, gradient=model$gr),
                      error = function(e) NA)
      
      # Record AIC of model fit
      aic.emp[i, disp.mod, mod] <- ifelse(!is.na(Opt[1]), 
                                          ifelse(sdreport(model)$pdHess, TMBhelper::TMBAIC(Opt),
                                          NA), NA)
      
      # Record model fit
      emp.mod.fits[[i]][[disp.mod]][[mod]] <- model
    }
    
  }
  
  # Setup data for "no dispersal" model
  to.fit <- make.inputs(dat, 'no.dispersal', 'neg.binom',
                        50, 'none')
  
  # Modify initial parameters for "no dispersal" model
  Parameters$log_sig_disp_mu <- 0
  Parameters$log_sig_disp_eps <- rep(0, dat$nperiods)
  
  # Fit "no dispersal" model
  model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_MM_sig", map=to.fit$Map,
                     random = to.fit$Random)
  model$env$beSilent()
  Opt <- tryCatch(nlminb(start=model$par, objective=model$fn, gradient=model$gr),
                  error = function(e) NA)
  
  # Record AIC of "no dispersal" model
  aic.no.dispersal[i] <- ifelse(!is.na(Opt[1]), 
                                    ifelse(sdreport(model)$pdHess, TMBhelper::TMBAIC(Opt),
                                           NA), NA)
  
  # Record "no dispersal" model fit
  no.dispersal.fits[[i]] <- model
  
}


# Plot Figure 2
################################################################################

# Build data
##################################

# Extract best model
# Only considering random effect and asymptote structures apparently
best.mods <- apply(aic.emp[,,3:4], 1, function(x) which(x==min(x, na.rm=TRUE), arr.ind=TRUE))
# best.mods[2,2] <- -1

# Setup empty container to hold predictions for each creek
report.ls <- list()

# Loop through creeks and extract predictions from best model
for(creek in 1:ncreeks) {
  report.ls[[creek]] <- emp.mod.fits[[creek]][[best.mods[1,creek]]][[best.mods[2,creek]+2]]$report()[[1]] %>% 
    data.frame()
  names(report.ls[[creek]]) <- c('obscount', 'predcount', 'times', 'distances', 
                                 'dist_factor')
  report.ls[[creek]]$creek <- paste('Creek', creek)
}

# Merge data between creeks
report <- do.call(rbind, report.ls)

# Summarize data
report_stats <- report %>%
  group_by(creek, times, distances) %>%
  summarize(obscount=mean(obscount), 
            predcount=first(predcount)) 

# Plot data
##################################

# Plot data
g1 <- ggplot(report_stats) +
  # Facet
  facet_wrap(~creek, scales = 'free') +
  # Plot data
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
             position=position_dodge(width=0.9), show.legend = FALSE) +
  # Labels
  labs(x='Time since release (days)', y='Number of recaptures') +
  # Legend
  scale_fill_discrete(name='Distance from\nrelease (m)') +
  # Theme
  ggsidekick::theme_sleek()
g1

# Export plot
ggsave(g1, filename=file.path(plotdir, "Fig2_obs_pred_data.png"), 
       width=7, height=5, units="in", dpi=600)


# Build Table 4
################################################################################

# Format results table ----------------------------------------------------

m.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='annual_mort')
t.max.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='t_max')
fifty.pct.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='fifty_pct_global')
fifty.pct.3wk.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='fifty_pct')
# Those indices should be unaffected by sampling structure for creek or 
# dispersal structure for model

table.width <- 4
print.mat <- matrix('', nrow=22, ncol=3 * table.width + 2)
print.mat[1,1] <- 'Random effect dispersal rate'
print.mat[1,table.width+2] <- 'Asymptotic dispersal rate'
print.mat[1,2*(table.width) + 3] <- 'Fixed effect dispersal rate'
row.marker <- 2

disp.arr <- mort.arr <- disp.3wk.arr <- array(NA, dim=c(3, 4, 3), 
                                              dimnames = list(disp.mod=c('normal', 'exponential', 'cauchy'),
                                                              creek=1:4,
                                                              mod.struc= c('fixed', 'random', 'asympBH')))

for(disp.mod in c('normal', 'exponential', 'cauchy')) {
  row.marker <- row.marker + 1
  print.mat[row.marker, c(1,table.width+2, 2*table.width+3)] <- Hmisc::capitalize(disp.mod)
  
  row.marker <- row.marker + 1
  print.mat[row.marker, c(1,table.width+2,2*table.width+3)] <- 'Location'
  print.mat[row.marker, c(2,table.width+3,2*table.width+4)] <- 'AIC'
  print.mat[row.marker, c(3,table.width+4,2*table.width+5)] <- 'M (1/yr)'
  print.mat[row.marker, c(4,table.width+5,2*table.width+6)] <- 'tmax (yrs)'
  print.mat[row.marker, c(5)] <- 'Mean dist @ 50% disp (m)'
  print.mat[row.marker, c(table.width+6)] <- 'Asymptotic dist @ 50% disp (m)'
  
  for(creek in 1:4) {
    row.marker <- row.marker + 1
    print.mat[row.marker, c(1, table.width+2, 2*table.width+3)] <- paste('Creek', creek)
    print.mat[row.marker, c(2, table.width+3, 2*table.width+4)] <- aic.emp[creek, disp.mod, 
                                                                           c('random effect', 
                                                                             'asymptote', 'fixed effect')] %>%
      round(1)
    for(mod.struc in 2:4) {
      if(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$pdHess) {
        
        print.mat[row.marker, c(NA,2*table.width+5, 3, table.width+4)[mod.struc]] <- paste0(
          round(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[m.ind], 2),
          '±',
          round(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$sd[m.ind], 2))
        mort.arr[disp.mod, creek, mod.struc-1] <- sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[m.ind]
        
        if(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[t.max.ind] < 100) {
          print.mat[row.marker, c(NA,2*table.width+6, 4, table.width+5)[mod.struc]] <- paste0(
            round(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[t.max.ind], 2),
            '±',
            round(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$sd[t.max.ind], 2))
        } else {
          print.mat[row.marker, c(NA,2*table.width+6, 4, table.width+5)[mod.struc]] <- paste(
            formatC(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[t.max.ind],
                    format='E', digits=2),
            '±',
            formatC(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$sd[t.max.ind], 
                    format='E', digits=2))
        }
        
        if(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[fifty.pct.ind] < 100) {
          print.mat[row.marker, c(NA,NA, 5, table.width+6)[mod.struc]] <- paste0(
            round(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[fifty.pct.ind], 2),
            '±',
            round(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$sd[fifty.pct.ind], 2))
        } else {
          print.mat[row.marker, c(NA,NA, 5, table.width+6)[mod.struc]] <- paste(
            formatC(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[fifty.pct.ind],
                    format='E', digits=2),
            '±',
            formatC(sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$sd[fifty.pct.ind], 
                    format='E', digits=2))
          
        }
        disp.arr[disp.mod, creek, mod.struc-1] <- sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[fifty.pct.ind]
        disp.3wk.arr[disp.mod, creek, mod.struc-1] <- sdreport(emp.mod.fits[[creek]][[disp.mod]][[mod.struc]])$value[fifty.pct.3wk.ind]
      }
    }
  }
  row.marker <- row.marker+1
}

# Remove results for models that did not converge
temp <- which(is.na(print.mat))#, arr.ind=TRUE)
print.mat[temp] <- 'Did not converge'
# print.mat[which(print.mat=='Did not converge')+1] <- ''

write.table(print.mat, file='original_model_cleaned/tables/model-fits.csv', row.names = FALSE, col.names = FALSE, sep=',')

mean(disp.arr[,1:2, 2:3], na.rm=T)
mean(disp.3wk.arr[,1:2, 2:3], na.rm=T)
mean(mort.arr[,1:2, 2:3], na.rm=T)

do.call(cbind, sdreport(emp.mod.fits[[2]]$exponential[[4]])[1:2])['sig_disp_sig',]
do.call(cbind, sdreport(emp.mod.fits[[2]]$exponential[[4]])[1:2])['sig_disp_sig',] * 0.95 / 0.05

