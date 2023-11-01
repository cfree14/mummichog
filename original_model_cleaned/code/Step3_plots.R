require(ggplot2)
require(gtable)
require(grid)
require(dplyr)
load('Code/sim-results.RData')

# Simulation estimates ----------------------------------------------------

# plot.sims <- function(rel.err.df, col, offset) {
#   dark.col <- stringr::str_replace(gplots::col2hex(col), 'FF', '99')
#   abline(v=offset, lty=3, col='black')
#   summary.stat.temp <- group_by(rel.err.df, est.mod) %>%
#     do(quants=quantile(.$rel.err, probs = c(.025, .25, .5, .75, .975), na.rm=TRUE))
#   summary.stat <- summary.stat.temp[match(c(disp.mods, 'no dispersal'), summary.stat.temp$est.mod),] %>%
#     filter(!is.na(est.mod))
#     
#   mapply(function(x,y) {
#     segments(x[1:2] + offset, y, x[5:4] + offset, y, lwd=2,
#              col=c('grey80', dark.col))
#     points(x[3] + offset, y, pch=21, col=dark.col, bg=col)}, 
#     summary.stat$quants, length(summary.stat$quants):1)
# }
# 
# max1=.3
# max2=.35
# max3=.75
# 
# # big.df <- res.df
# # res.df <- filter(res.df, est.mod != 'no dispersal')
# 
# png('Figs/poisson-sims.png', height=5.5, width=4.8, units = 'in', res=200)
# par(mfrow=c(n_distinct(res.df$sim.mod),1), mar=rep(.5, 4), oma=c(5,6,4,1))
# plot.max <- n_distinct(res.df$est.mod) + 0.5
# for(mod in disp.mods) {
#   plot(1,1, xlim=c(-max1, max1+max2+2*max3), ylim=c(.5, plot.max), 
#        type='n', axes=FALSE, ann=FALSE)
#   text(-.7, plot.max+0.2, Hmisc::capitalize(mod), xpd=NA, pos=4)
#   text(-.63, (plot.max-0.5):1, Hmisc::capitalize(c(disp.mods, 'no dispersal')), xpd=NA, pos=4)
#   if(mod==disp.mods[1]) {
#     text(c(0, max1+max2, max1+max2+max3), plot.max, xpd=NA, pos=3, offset=2,
#          c('T max', '% w/in 50m', 'Dist @ 50%'),
#          col=stringr::str_replace(gplots::col2hex(c('red', 'blue', 'green')), 'FF', '99'))
#   }
#   
#   temp <- filter(res.df, sim.mod==mod)  %>%
#     mutate(rel.err = (t_max-t.max) / t.max) %>%
#     plot.sims(col='red', offset=0)
#     
#   temp <- filter(res.df, sim.mod==mod)  %>%
#     mutate(rel.err = (pct_at_dist-true.val.mat[sim.mod,'pct_at_dist']) /
#              true.val.mat[sim.mod,'pct_at_dist']) %>%
#    plot.sims(col='blue', offset=max1+max2)
# 
#   temp <- filter(res.df, sim.mod==mod)  %>%
#     mutate(rel.err = (fifty_pct-true.val.mat[sim.mod,'fifty_pct']) /
#              true.val.mat[sim.mod,'fifty_pct']) %>%
#     plot.sims(col='green', offset=max1+max2+max3)
# }
# box('inner')
# axis(1, at=max1+max2+max3+c(0, .25, .5), labels = c('', 0.25, 0.5), line = 1)
# axis(1, at=max1+max2+c(0, .2, .4), labels=c('', .2, .4), line=1)
# axis(1, at=c(0, .2), labels=c('', .2), line=1)
# mtext('Relative error', 1, 3, outer=TRUE)
# dev.off()
# 
# system('open Figs/poisson-sims.png')


# Survival estimates ----------------------------------------------------
annual.mort.true <- annual.mort

sim.res.df$est.mod[sim.res.df$est.mod=='no.dispersal'] <- ''                                                                       
sim.res.df$disp.str[sim.res.df$disp.str=='asympBH'] <- 'asympMM'                                                                       

this.plot <- sim.res.df %>%
  mutate(rel.err = (annual.mort-annual.mort.true)/annual.mort.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05, na.rm=T), low=quantile(rel.err, .25, na.rm=T), 
            mid=median(rel.err, na.rm=T), high=quantile(rel.err, .75, na.rm=T), 
            # max=ifelse(quantile(rel.err, .95)<.5, quantile(rel.err, .95), .5)) %>%
            max=quantile(rel.err, .95, na.rm=T)) %>%
  ggplot() +
  geom_point(aes(x=mid, y=disp.str), cex=2) +
  geom_segment(aes(x=min, xend=max, y=disp.str, yend=disp.str)) +
  geom_segment(aes(x=low, xend=high, y=disp.str, yend=disp.str), lwd=1.5) +
  geom_vline(xintercept=0, lty=2) +
  facet_grid(est.mod ~ sim.mod, scales = 'free_y', space='free_y') +
  ggsidekick::theme_sleek() +
  theme(strip.text.y = element_text(color = 'grey30'), 
        strip.text.x = element_text(color = 'grey30')) +
  xlab('Relative error in M') + ylab('Dispersal structure')

z <- ggplotGrob(this.plot)
posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z$layout, grepl("strip-t", name), select = t:r)

width <- z$widths[max(posR$r)]    # width of current right strips
height <- z$heights[min(posT$t)]  # height of current top strips

z <- gtable_add_cols(z, width, max(posR$r))  %>% 
  gtable_add_rows(height, min(posT$t)-1)

stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, fill = NA)),
  textGrob('Estimation dispersal distribution', rot = -90, gp = gpar(fontsize = 11, col = "grey30"))))
stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, fill = NA)),
  textGrob('Simulation dispersal distribution', gp = gpar(fontsize = 11, col = "grey30"))))

z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right") %>%
  gtable_add_grob(stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

png('original_model_cleaned/figures/survival-results.png', height=8, width=8, units = 'in', res=200)
grid.newpage()
grid.draw(z)
dev.off()
system('open original_model_cleaned/figures/survival-results.png')

sim.res.df %>%
  mutate(rel.err = (annual.mort-annual.mort.true)/annual.mort.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05, na.rm=T), low=quantile(rel.err, .25, na.rm=T), 
            mid=median(rel.err, na.rm=T), high=quantile(rel.err, .75, na.rm=T), 
            max=quantile(rel.err, .95, na.rm=T)) %>%
  # filter(est.mod != 'no dispersal', disp.str != 'constant') %>%
  arrange(sim.mod, desc(abs(mid))) %>% 
  filter(est.mod != '', disp.str != 'no dispersal') %>% 
  arrange(desc(abs(mid)))
  group_by(disp.str) %>% summarize(mean(mid, na.rm=T))

# dispersal estimates ----------------------------------------------------
sim.res.df$fifty.pct.true <- sapply(sim.res.df$sim.mod, function(x)
  switch(x,
         normal = qnorm(.75, 0, sig.disp[4]),
         exponential = qexp(.5, 1/sig.disp[4]),
         cauchy = qcauchy(.75, sig.disp[4])))

this.plot <- sim.res.df %>%
  filter(est.mod != 'no dispersal', est.mod != 'no.dispersal', est.mod != '') %>%
  mutate(rel.err = (fifty.pct - fifty.pct.true)/fifty.pct.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05, na.rm=T), low=quantile(rel.err, .25, na.rm=T), 
            mid=median(rel.err, na.rm=T), high=quantile(rel.err, .75, na.rm=T), 
            max=quantile(rel.err, .95, na.rm=T)) %>% #ifelse(quantile(rel.err, .95)<.5, quantile(rel.err, .95), .5)) %>%
  ggplot() +
  geom_point(aes(x=mid, y=disp.str), cex=2) +
  geom_segment(aes(x=min, xend=max, y=disp.str, yend=disp.str)) +
  geom_segment(aes(x=low, xend=high, y=disp.str, yend=disp.str), lwd=1.5) +
  geom_vline(xintercept=0, lty=2) +
  facet_grid(est.mod ~ sim.mod, scales = 'free_y', space='free_y') +
  ggsidekick::theme_sleek() +
  xlab('Relative error in distance at 50% dispersal on day 20') +
  ylab('Dispersal structure')

z <- ggplotGrob(this.plot)
posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z$layout, grepl("strip-t", name), select = t:r)

width <- z$widths[max(posR$r)]    # width of current right strips
height <- z$heights[min(posT$t)]  # height of current top strips

z <- gtable_add_cols(z, width, max(posR$r))  %>% 
  gtable_add_rows(height, min(posT$t)-1)

stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, fill = NA)),
  textGrob('Estimation dispersal distribution', rot = -90, gp = gpar(fontsize = 11, col = "grey30"))))
stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, fill = NA)),
  textGrob('Simulation dispersal distribution', gp = gpar(fontsize = 11, col = "grey30"))))

z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right") %>%
  gtable_add_grob(stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

png('original_model_cleaned/figures/dispersal-results.png', height=8, width=8, units = 'in', res=200)
grid.newpage()
grid.draw(z)
dev.off()
system('open original_model_cleaned/figures/dispersal-results.png')

sim.res.df %>%
  mutate(rel.err = (fifty.pct-fifty.pct.true)/fifty.pct.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05, na.rm=T), low=quantile(rel.err, .25, na.rm=T), 
            mid=median(rel.err, na.rm=T), high=quantile(rel.err, .75, na.rm=T), 
            max=quantile(rel.err, .95, na.rm=T)) %>%
  # filter(est.mod != 'no dispersal', disp.str != 'constant') %>%
  arrange(sim.mod, desc(abs(mid))) %>% 
  # filter(est.mod != '', disp.str != 'no dispersal') %>% 
  # arrange(desc(abs(mid)))
group_by(disp.str) %>% summarize(mean(mid, na.rm=T))


# Simulation AIC ----------------------------------------------------------

aic.winner <- apply(aic.sim[,1:3,], c(1,3), which.min) 
winning.pct <- numeric(3)
for(ii in 1:3) {
  winning.pct[ii] <- sum(aic.winner[ii,] == ii)/nreps
}

aic.mns <- apply(aic.sim[,1:3,], 1:2, mean, na.rm=TRUE) 
aic.sds <- apply(aic.sim[,1:3,], 1:2, sd)
disp.mods <- rownames(aic.mns)

paste(round(aic.mns, 1), '±', round(aic.sds, 1)) %>% c(paste0(winning.pct*100, '%')) %>% 
  matrix(nrow=3, dimnames=list('Simulation model'=Hmisc::capitalize(disp.mods), 
                               c(Hmisc::capitalize(disp.mods), 
                                 '% correct model chosen'))) %>%
  write.csv(file='original_model_cleaned/figures/AIC.csv')


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

write.table(print.mat, file='original_model/Figs/model-fits.csv', row.names = FALSE, col.names = FALSE,
            sep=',')

mean(disp.arr[,1:2, 2:3], na.rm=T)
mean(disp.3wk.arr[,1:2, 2:3], na.rm=T)
mean(mort.arr[,1:2, 2:3], na.rm=T)

do.call(cbind, sdreport(emp.mod.fits[[2]]$exponential[[4]])[1:2])['sig_disp_sig',]
do.call(cbind, sdreport(emp.mod.fits[[2]]$exponential[[4]])[1:2])['sig_disp_sig',] * 0.95 / 0.05


# Plot data and fits from best model --------------------------------------

require(ggplot2)

best.mods <- apply(aic.emp[,,3:4], 1, function(x) which(x==min(x, na.rm=TRUE), arr.ind=TRUE))
# best.mods[2,2] <- -1
report.ls <- list()

for(creek in 1:4) {
  report.ls[[creek]] <- emp.mod.fits[[creek]][[best.mods[1,creek]]][[best.mods[2,creek]+2]]$report()[[1]] %>% 
    data.frame()
  names(report.ls[[creek]]) <- c('obscount', 'predcount', 'times', 'distances', 
                                 'dist_factor')
  report.ls[[creek]]$creek <- paste('Creek', creek)
}

report <- do.call(rbind, report.ls)

png('original_model/Figs/obs-pred-data.png', height=5, width=7, units = 'in', res=200)
report %>%
  group_by(creek, times, distances) %>%
  summarize(obscount=mean(obscount), predcount=first(predcount)) %>%
  ggplot() +
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
             position=position_dodge(width=0.9), show.legend = FALSE) +
  facet_wrap(~creek, scales = 'free') +
  xlab('Time since release (days)') +
  ylab('Number of recaptures') +
  scale_fill_discrete(name='Distance from\nrelease (m)') +
  ggsidekick::theme_sleek()
dev.off()
system('open original_model/Figs/obs-pred-data.png')

pdf('original_model/Figs/obs-pred-data-talk.pdf')

ggplot() + theme_classic()

report %>%
  filter(creek=='Creek 2', times==2) %>%
  # group_by(times, distances) %>%
  group_by(distances) %>%
  summarize(obscount=mean(obscount), predcount=first(predcount), 
            times=factor(2, levels = unique(report$times[report$creek=='Creek 2']))) %>%
  ggplot() +
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  # geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
  #            position=position_dodge(width=0.9), show.legend = FALSE) +
  xlab('Time since release (days)') +
  ylab('Number of recaptures') +
  scale_fill_discrete(name='Distance from\nrelease (m)', drop=FALSE) +
  scale_x_discrete(limits=factor(unique(report$times[report$creek=='Creek 2'])), drop=FALSE) +
  ylim(obscount=1.05*c(0,max(report$obscount))) +
  theme_classic(base_size = 18)

report %>%
  filter(creek=='Creek 2') %>%
  group_by(times, distances) %>%
  # group_by(distances) %>%
  summarize(obscount=mean(obscount), predcount=first(predcount)) %>%
  ggplot() +
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  # geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
  #            position=position_dodge(width=0.9), show.legend = FALSE) +
  xlab('Time since release (days)') +
  ylab('Number of recaptures') +
  scale_fill_discrete(name='Distance from\nrelease (m)', drop=FALSE) +
  ylim(obscount=1.05*c(0,max(report$obscount))) +
  theme_classic(base_size = 18)
dev.off()


# Summary data for appendix table -----------------------------------------

report.new <- mutate(report, creek = as.numeric(gsub('Creek ', '', creek)))
library(lubridate)

dates <- c('10/11/12', '10/12/12', '10/13/12', '10/14/12', '10/20/12', '11/1/12', '12/8/12', '5/20/12', '5/21/12', 
           '3/12/13', '6/2/13', '7/3/13', '7/27/13', '9/14/13') %>% mdy

tagging.dates <- dates[c(2,1,8,9)]
out.list <- list()

for(cr in 1:4) {
  dat <- filter(report.new, creek == cr)
  times <- unique(dat$times)
  out.list[[cr]] <- tibble(date = tagging.dates[cr] + ddays(times)) %>%
    cbind(times) %>%
    right_join(dat) %>%
    group_by(creek, times) %>%
    summarize(count = sum(obscount))
}
