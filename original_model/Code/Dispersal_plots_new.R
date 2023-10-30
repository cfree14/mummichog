##plotting output from the dispersal model - Olaf Jensen Aug 25, 2015

#read in the output matrix
data=read.table("C:/Data/Rutgers/Gomex_fundulus/Dispersal/ADMB/dm_normal_dyn_sig.rep", header=TRUE)

#plot overall observations vs. predictions 
plot(data$predcount, data$obscount, ylab = "Observed count", xlab = "Predicted Count")

#subset data
unique_times=unique(data$times)

for(i in 1:9)
{
  data2=subset(data,times==unique_times[i])
  
  #plot observations vs. distance
  plot(data2$distances,data2$obscount, ylim=range(c(0,1.1*max(data2$obscount))),ylab = "count", xlab = "distance (m)", main=unique_times[i])
  #plot(data2$distances,data2$obscount, ylim=range(c(0,1.1*max(data2$obscount))),ylab = "count", xlab = "distance (m)")
  
  
  #plot predictions vs. distance on the same axes
  par(new = TRUE)
  plot(data2$distances, data2$predcount, ylim=range(c(0,1.1*max(data2$obscount))), axes = FALSE, xlab = "", ylab = "", type="l")
}

##Plotting change in dispersal variance through time
