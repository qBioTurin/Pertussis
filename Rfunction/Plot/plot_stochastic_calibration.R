library(ggplot2)
library(parallel)
library(cowplot)

par_processors <- 16
# Plot the optimal trace
n_sim <- 2^12
optim_traces <- read.csv("~/results_stochastic_model_analysis/Pertussis-analysys-1.trace", header = TRUE, sep = "")

hyst <- read.csv("~/input/reference_data.csv", header = FALSE, sep = "")
hyst <- data.frame(Trace=unlist(hyst), Time=c(0:(length(unlist(hyst))-1))+1974)

y_names <-names(optim_traces)
# Get the indexes of all infected places
indexes.p <- c( which(y_names %in% grep("Ip_*", y_names, value=T)) )
indexes.s <- c( which(y_names %in% grep("Is_*", y_names, value=T)) )
indexes<-c(indexes.p,indexes.s)
indexes.counters<-c( which(y_names %in% grep("InfectCount_*", y_names, value=T)) )
# Get the time interval from the simulations
time<-unique(optim_traces$Time)
# Compute the increments in the number of infects year by year
infect<-sapply(0:(n_sim-1),
               function(k){
                 out<-optim_traces[(2+k*(length(time))):(length(time)+k*length(time)),indexes.counters]
                 infected<-rowSums(out)
                 c(infected[1],diff(infected,differences = 1))
               }
)

# Give an unique id to each simulation
sim_optim_traces<-data.frame(Time=rep( time[-c(1)]/365+1973,n_sim),infected=c(infect),ID=rep(1:n_sim,each=(length(time)-1) ))

# PLOT 1: all traces of the best fitting experiment, their mean value and the hystorical data

# Comupte mean, median and the standard deviation
mean<-apply(infect,1,mean)
median <- apply(infect,1,median)
sd<-apply(infect,1,sd)
# Compute the "confidence interval" for the mean
up<-mean+(qnorm(.975))*sd/sqrt(n_sim)
lp<-mean-(qnorm(.975))*sd/sqrt(n_sim)

stats<-data.frame(days = time[-c(1)]/365+1973,mean=mean,up=up,lp=lp)


g.boxplot<-ggplot()+
  geom_boxplot(data = sim_optim_traces,aes(x=Time,y=infected,group=Time))+
  geom_line(data = hyst, aes(x=Time,y=Trace, col="RealInfects")) +
  geom_line(data = stats, aes(x=days,y=mean, col="Mean"))+
  labs(x="Years",y="Number of cases")+
  scale_color_manual("",
                     breaks=c("RealInfects","Mean"),
                     labels=c("Real cases","Mean"),
                     values = c("red","blue"),
                     limits=c("RealInfects","Mean"))+
  theme(
    axis.text=element_text(size=18),
    axis.title=element_text(size=20,face="bold"),
    legend.text=element_text(size=18),
    legend.title=element_text(size=20,face="bold"),
    legend.position="bottom",
    legend.box = "vertical",
    legend.key.size = unit(1.3, "cm"),
    legend.key.width = unit(1.3,"cm"))+
    scale_x_continuous(limits=c(1973.5,1994.5),breaks = seq(1974,1994,5) )

# Save figures to file ( boxplot and trajectories for all the traces)
ggsave(g.boxplot, file="~/Plots/all_traces_Boxplot.pdf", device="pdf", dpi="retina", width=8, height=6, units="in")

#
folder <- "~/results_stochastic_model_calibration/"

n_sim <- 256

rank <- read.csv(paste0(folder,"Pertussis-calibration_optim-trace.csv"), sep = "")
rank <- rank[order(rank$distance),]

traces <- lapply(rank$id[1:9],function(x){
  trace <- read.csv(paste0(folder,"Pertussis-calibration-",x,".trace"), header = TRUE, sep = "")
  
  # Compute the increments in the number of infects year by year
  infect<-sapply(0:(n_sim-1),
                 function(k){
                   out<-trace[(2+k*(length(time))):(length(time)+k*length(time)),indexes.counters]
                   infected<-rowSums(out)
                   c(infected[1],diff(infected,differences = 1))
                 }
  )
  sim_trace<-data.frame(Time=rep( time[-c(1)]/365+1973,n_sim),infected=c(infect),EXP=x,ID=rep(1:n_sim,each=(length(time)-1) ))
  return(sim_trace)
})

# Add the best fitting experiment
traces <- do.call("rbind",traces)
op_tr <- sim_optim_traces[1:n_sim,]
op_tr$EXP <- rep(min(traces$EXP)-1,n_sim)
traces <- rbind(traces, op_tr)

cl <- makeForkCluster(nnodes = par_processors)
# Compute year by year statistics for each experiment  
stats<-parLapply(X=unique(traces$EXP),
                 function(k,tr){
                            infect<-tr[which(tr$EXP==k), ] #-22]
                            infect_day <- split(infect$infected,infect$Time)
                            mean <- as.data.frame(sapply(infect_day,mean))
                            median <- as.data.frame(sapply(infect_day,median))
                            sd <- as.data.frame(sapply(infect_day,sd))
                            # Compute the "confidence interval" for the mean
                            up<-mean+(qnorm(.975))*sd/sqrt(n_sim)
                            lp<-mean-(qnorm(.975))*sd/sqrt(n_sim)
                            data.frame(days = time[-c(1)]/365+1973,mean=mean,up=up,lp=lp,IDconf=paste(k))
                          },
                 cl = cl, 
                 tr = traces)

stats<-do.call("rbind", stats)
names(stats) <- c("days","mean","up","lp","IDconf")

Meanarea <- parLapply(X=unique(stats$days),
                      function(x, stats)
                      {
                        minMean<- min(stats$mean[which(stats$days==x)])
                        maxMean<-max(stats$mean[which(stats$days==x)])
                        minSD<- min(stats$lp[which(stats$days==x)])
                        maxSD<-max(stats$up[which(stats$days==x)])
                        c(x,minMean,maxMean,minSD,maxSD)
                      },
                      stats=stats,
                      cl=cl)

stopCluster(cl = cl)

Meanarea<-as.data.frame(do.call("rbind", Meanarea))

colnames(Meanarea)<-c("time","min","max","lp","up")


g <- ggplot()+
      geom_line(data = traces,aes(x=Time,y=infected,group=ID),col="grey")+
      geom_ribbon(data=Meanarea, 
              aes(x=time,ymin=lp,ymax=up), fill="green",col="green", alpha=.1,linetype="dashed")+
      geom_ribbon(data=Meanarea, 
              aes(x=time,ymin=min,ymax=max), fill="blue", alpha=.6)+
      geom_line(data = hyst, aes(x=Time,y=Trace, col="Real cases")) +
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"),
            legend.text=element_text(size=18),
            legend.title=element_text(size=20,face="bold"),
            legend.position="bottom",
            legend.key.size = unit(1.3, "cm"),
            legend.key.width = unit(1.3,"cm") )+
    scale_color_manual("",breaks=c("Real cases","Mean","Simulations","Standard Deviation"),labels=c("Real cases","Mean","Simulations","Standard Deviation"),values = c("red","blue","grey","green"),limits=c("Real cases","Mean","Simulations","Standard Deviation"))+
    labs(x="Years", y="Number of cases")+
    scale_x_continuous(limits=c(1974,1994),breaks = seq(1974,1994,5) )

ggsave(g, file="~/Plots/top_n_traces.pdf", device="pdf", dpi="retina", width=8, height=6, units="in")


plt <- plot_grid(g,g.boxplot,
                 labels = c("a)", "b)"),
                 ncol = 1, nrow = 2)

ggsave(plt, file="~/Plots/CalibrationStoch.pdf", device="pdf", dpi="retina", width=10, height=15, units="in")
