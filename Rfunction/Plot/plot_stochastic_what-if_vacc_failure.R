library(ggplot2)
library(cowplot)
#######################
# Configuration:
experiment <- 2
folder <- paste0("~/results_model_analysis_vacc_failure_0",experiment,"/")
fname <- "Pertussis-analysys-1.trace"
surveillance <- "~/input/reference_data_ext.csv"
#######################

filter<-function(output)
{
  ynames <- names(output)
  col_names <- "(InfectCount_a){1}[0-9]{1}"
  col_idxs <- c( which(ynames %in% grep(col_names, ynames, value=T)) )
  # Reshape the vector to a row vector
  time <- unique(output$Time)
  ret <- lapply(0:(length(output[,1])/length(time)-1),
                  function(k){
                    out<-output[(2+k*(length(time))):(length(time)+k*length(time)),col_idxs]
                    infected<-rowSums(out)
                    data.frame(Time = time[-1],
                               infects = c(infected[1],diff(infected,differences = 1)),
                               id=rep(k+1,length(time)-1))
                  }
  )
  ret <- do.call("rbind",ret)
  return(ret)
}

trace <- filter(read.csv(paste0(folder,fname), sep = "", header=TRUE))
hyst <- as.data.frame(t(read.csv(surveillance, header = FALSE, sep = "")))
hyst$Time <- unique(trace$Time)
names(hyst) <- c("Surveillance","Time")
mean <- as.data.frame(aggregate(trace$infects, list(trace$Time), mean))
names(mean) <- c("Time","Mean")
g <- ggplot() +
  geom_line(data=trace, aes(y=infects, x=Time/365+1973, group=id), col="grey") +
  geom_line(data=mean, aes(y=Mean, x=Time/365+1973), col="blue") +
  geom_line(data=hyst, aes(y=Surveillance, x=Time/365+1973), col ="red", size=1.05)+
  labs(x="Years",y="Number of cases")+
  scale_color_manual("",
                     breaks=c("RealInfects","Mean","Simulations"),
                     labels=c("Real cases","Mean","Simulations"),
                     values = c("red","blue","grey"),
                     limits=c("RealInfects","Mean","Simulations"))+
  theme(
    axis.text=element_text(size=18),
    axis.title=element_text(size=20,face="bold"),
    legend.text=element_text(size=18),
    legend.title=element_text(size=20,face="bold"),
    legend.position="bottom",
    legend.box = "vertical",
    legend.key.size = unit(1.3, "cm"),
    legend.key.width = unit(1.3,"cm"))+
  scale_x_continuous(limits=c(1973.5,2016.5),breaks = seq(1974,2016,5) )

g.boxplot <- ggplot() +
  geom_boxplot(data=trace, aes(y=infects, x=Time/365+1973, group=Time)) +
  geom_line(data=mean, aes(y=Mean, x=Time/365+1973), col="blue") +
  geom_line(data=hyst, aes(y=Surveillance, x=Time/365+1973), col ="red", size=1.05)+
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
  scale_x_continuous(limits=c(1973.5,2016.5),breaks = seq(1974,2016,5) )

g.boxplot.zoom <- ggplot() +
  geom_boxplot(data=trace, aes(y=infects, x=Time/365+1973, group=Time)) +
  geom_line(data=mean, aes(y=Mean, x=Time/365+1973), col="blue") +
  geom_line(data=hyst, aes(y=Surveillance, x=Time/365+1973), col ="red", size=1.05)+
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
  scale_x_continuous(limits=c(1993.5,2016.5),breaks = seq(1994,2016,5) )

plt <- plot_grid(g,g.boxplot,g.boxplot.zoom,
                 labels = c("a)", "b)", "c)"),
                 ncol = 1, nrow = 3)

ggsave(plt, file=paste0("~/Plots/prob0",experiment,".pdf"), device="pdf", dpi="retina", width=10, height=22.5, units="in")