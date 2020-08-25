
library(ggplot2)
library(viridis)

rank <- read.csv("Results/results_model_calibration/Pertussis-calibration_optim-trace.csv", header = TRUE, sep = "")


names<-c("prob_boost","prob_infectionS","prob_infectionR_l1",
         "init_S_a1","init_S_a2","init_S_a3","init_R_a1_nv_l4",
         "init_R_a2_nv_l1", "init_R_a2_nv_l2", "init_R_a2_nv_l3", "init_R_a2_nv_l4",
         "init_R_a3_nv_l1", "init_R_a3_nv_l2", "init_R_a3_nv_l3", "init_R_a3_nv_l4")

rank <- rank[order(rank$distance,decreasing = F),]

optim<-rank[1,-c(1:2)]
colnames(optim)<-names
bestconfig<-optim

cat("Optimal configuration:\n")
optim[1,]

  
  filter<-function(output)
  {
    ynames <- names(output)
    col_names <- "(InfectCount_a){1}[0-9]{1}"
    col_idxs <- c( which(ynames %in% grep(col_names, ynames, value=T)) )
    # Reshape the vector to a row vector
    ret <- rowSums(output[,col_idxs])[-1]
    ret <- c(ret[1],diff(ret, differences = 1))
    return(ret)
  }
  
  traces <- lapply(rank$id,function(x){
    if(file.exists(paste0("Results/results_model_calibration/Pertussis-calibration-",x,".trace"))){
      trace <- unlist(filter(read.csv(paste0("Results/results_model_calibration/Pertussis-calibration-",x,".trace"), sep = "")))
      ids <- rep(x, length(trace))
      rnk <- rep((rank$distance[which(rank$id==x)]-min(rank$distance))/max(rank$distance), length(trace))
      trace <- data.frame(Trace=trace,id=ids, Time=1974+c(0:(length(unlist(trace))-1)), Score=1-rnk )
      return(trace)
    }
  })
  
  traces <- as.data.frame(do.call("rbind",traces))
  
  
  trace <- filter(read.csv(paste0("Results/results_model_calibration/Pertussis-calibration-",rank[1,]$id,".trace"), sep = ""))
  
  trace <- data.frame(Trace=trace, Time=1974+c(0:(length(unlist(trace))-1)))
  
  hyst <- read.csv("./input/reference_data_ext.csv", header = FALSE, sep = "")[1:21]
  
  hyst <- data.frame(Trace=t(as.data.frame(hyst)), Time=1974+c(0:(length(unlist(hyst))-1)))
  traces <- traces[order(traces$Score, decreasing = FALSE),]
  
  indexesToHigh<-traces[which(traces$Trace>78000),]$id
  
  summary(traces$id)->s
  moreindexes<-c(sample(1:s[5],65000),sample(s[5]:s[6],5000) )
  
  indexesTodelete<-which(traces$id %in% unique(c(indexesToHigh,moreindexes)) )
  
  plt<-ggplot()+
    geom_line(data=traces[-indexesTodelete,],aes(y=Trace,group=id,col=1-(Score-min(Score))/(max(Score)-min(Score)),x=Time)) +
    geom_line(data=trace,aes(y=Trace,x=Time),col="black", size = 1.2, linetype = 1) +
    geom_line(data=hyst,aes(y=Trace,x=Time),col="orangered2", size = 1.2, linetype = 1) +
    theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold"))+
    labs(x="Years", y="Number of cases")+
    scale_colour_gradientn("Error",
                           colours = c("#FFFF00","#33CC33","#330066"),
                           values= c(0,.07,1),
                           breaks=c(0,1),
                           labels=c("min","max"))+
    scale_x_continuous(limits=c(1974,1994),breaks = seq(1974,1994,5) )
  