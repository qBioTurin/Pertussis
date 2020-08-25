library(ggplot2)
library(cowplot)


rank <- read.csv("Results/results_model_calibration_STOCHASTIC/Pertussis-calibration_optim-trace.csv", header = TRUE, sep = "")
rank <- rank[order(rank$distance),]
rank$id[1:10]->id

names<-c("prob_boost","prob_infectionS","prob_infectionR_l1",
         "init_S_a1","init_S_a2","init_S_a3","init_R_a1_nv_l4",
         "init_R_a2_nv_l1", "init_R_a2_nv_l2", "init_R_a2_nv_l3", "init_R_a2_nv_l4",
         "init_R_a3_nv_l1", "init_R_a3_nv_l2", "init_R_a3_nv_l3", "init_R_a3_nv_l4")


optim<-rank[1,-c(1:2)]
colnames(optim)<-names
paste("Optimal configuration=")
optim.stoch=optim[1,]
optim.stoch

  
  rank <- read.csv("Results/results_model_calibration_STOCHASTIC/Pertussis-calibration_optim-trace.csv", header = TRUE, sep = "")
  rank <- rank[order(rank$distance),]
  rank$id[1:10]->id
  
  list.files("./Results/results_model_calibration_STOCHASTIC/", pattern = ".trace")-> listFile
  
  ind<-which(listFile %in% paste0("Results/results_model_calibration_STOCHASTIC/Pertussis-calibration-",id,".trace")) 
  
  n_sim <- 40
  length(listFile[ind])->n_conf
  
  trace<-lapply(listFile[ind],
                function(x){
                  return(read.csv(x, sep = ""))
                })
  length(trace[[1]][,1])->l
  trace <- do.call("rbind", trace)
  
  hyst <- read.csv("./input/reference_data_ext.csv", header = FALSE, sep = "")[1:21]
  hyst <- data.frame(Trace=unlist(hyst), Time=c(0:(length(unlist(hyst))-1))+1974)
  
  y_names <-names(trace)
  # Get the indexes of all infected places
  indexes.p <- c( which(y_names %in% grep("Ip_*", y_names, value=T)) )
  indexes.s <- c( which(y_names %in% grep("Is_*", y_names, value=T)) )
  indexes<-c(indexes.p,indexes.s)
  indexes.counters<-c( which(y_names %in% grep("InfectCount_*", y_names, value=T)) )
  # Get the time interval from the simulations
  time<-unique(trace[,1])
  # Compute the increments in the number of infects year by year
  infect<-lapply(0:((n_sim)*n_conf-1),
                 function(k){
                   
                   out<-trace[(2+k*(length(time))):(length(time)+k*length(time)),indexes.counters]
                   infected<-rowSums(out)
                   c(infected[1],diff(infected,differences = 1),id[as.integer(k/n_sim)+1])
                   
                 }
  )
  
  infect2 <- do.call("rbind", infect)
  colnames(infect2) <- c(paste0("time",1:21),"ID")
  
  infect <-c(infect2[,-22])
  
  # Give an unique id to each simulation
  #sim_traces<-data.frame(Time=rep( time[-c(1)]/365+1973,n_sim*n_conf),infected=c(infect),ID=rep(1:n_sim,each=(length(time)-1) ))
  
  sim_traces<-data.frame(Time=rep( time[-c(1)]/365+1973,n_sim*n_conf ),infected=c(infect),IDconf=rep(id,each= 21*n_sim ),IDsim=rep(1:(n_conf*n_sim),each= length(time[-1])) )
  
  # PLOT 1: all traces of the same experiment, their mean value and the hystorical data
  # Comupte mean, median and the standard deviation
  stats<-lapply(id,function(k)
  {
    infect3<-infect2[which(infect2[,"ID"]==k),-22]
    mean<-apply(infect3,2,mean)
    median <- apply(infect3,2,median)
    sd<-apply(infect3,2,sd)
    # Compute the "confidence interval" for the mean
    up<-mean+(qnorm(.975))*sd/sqrt(n_sim)
    lp<-mean-(qnorm(.975))*sd/sqrt(n_sim)
    
    data.frame(days = time[-c(1)]/365+1973,mean=mean,up=up,lp=lp,IDconf=paste(k))
  }
  )
  
  stats<-do.call("rbind", stats)
  
  Meanarea<-t(sapply(unique(stats$days), function(x){
    minMean<- min(stats$mean[which(stats$days==x)])
    maxMean<-max(stats$mean[which(stats$days==x)])
    minSD<- min(stats$lp[which(stats$days==x)])
    maxSD<-max(stats$up[which(stats$days==x)])
    c(x,minMean,maxMean,minSD,maxSD)
  }))
  Meanarea<-as.data.frame(Meanarea)
  colnames(Meanarea)<-c("time","min","max","lp","up")
  
  summary(rank[1:10,-(1:2)])->summarySimul
  summarySimul[c(1,6),]
  
  
  pl1<-ggplot()+
    geom_line(data = sim_traces,aes(x=Time,y=infected,group=IDsim),col="grey")+
    geom_line(data = hyst, aes(x=Time,y=Trace, col="Real cases")) +
    #geom_line(data = stats, aes(x=days, y=mean, col=IDconf, group=IDconf ))+
    geom_ribbon(data=Meanarea, 
                aes(x=time,ymin=lp,ymax=up), fill="green",col="green", alpha=.1,linetype="dashed")+
    geom_ribbon(data=Meanarea, 
                aes(x=time,ymin=min,ymax=max), fill="blue", alpha=.6)+
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
  
  
  
  pl2<-ggplot()+
    geom_boxplot(data = sim_traces[which(sim_traces$IDconf==rank[1,]$id),],aes(x=Time,y=infected,group=Time))+
    geom_line(data = hyst, aes(x=Time,y=Trace, col="Real cases")) +
    geom_line(data = stats[which(stats$IDconf==rank[1,]$id),], aes(x=days,y=mean, col="Mean"))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=18),
          legend.title=element_text(size=20,face="bold"),
          legend.position="bottom",
          legend.key.size = unit(1.3, "cm"),
          legend.key.width = unit(1.3,"cm") )+
    scale_color_manual("",breaks=c("Real cases","Mean"),labels=c("Real cases","Mean"),values = c("red","blue"),limits=c("Real cases","Mean"))+
    labs(x="Years", y="Number of cases")+
    scale_x_continuous(limits=c(1973,1995),breaks = seq(1974,1994,5) )
  
  plot_grid(pl1,pl2,
            labels = c("a)", "b)"),
            ncol = 1, nrow = 2)