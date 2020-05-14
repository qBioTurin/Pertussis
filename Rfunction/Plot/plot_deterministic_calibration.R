library(ggplot2)
library(viridis)
library(parallel)

folder <- "~/results_deterministic_model_calibration/"

par_processors <- 16

rank <- read.csv(paste0(folder,"Pertussis-calibration_optim-trace.csv"), header = TRUE, sep = "")

names<-c("prob_boost","prob_infectionS","prob_infectionR_l1",
         "init_S_a1","init_S_a2","init_S_a3","init_R_a1_nv_l4",
         "init_R_a2_nv_l1", "init_R_a2_nv_l2", "init_R_a2_nv_l3", "init_R_a2_nv_l4",
         "init_R_a3_nv_l1", "init_R_a3_nv_l2", "init_R_a3_nv_l3", "init_R_a3_nv_l4")


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

optim<-rank[43044,-c(1:2)]
rank <- head(rank,43044)
rank <- rank[order(rank$distance,decreasing = F),]


colnames(optim)<-names
bestconfig<-optim

cat("Optimal configuration:\n")
optim[1,]

cl <- makeForkCluster(nnodes = par_processors)

traces <- parLapply(X=rank$id,
                    function(x,rnk,dir)
                      {
                          if(file.exists(paste0(dir,"Pertussis-calibration-",x,".trace")))
                            {
                              trace <- unlist(filter(read.csv(paste0(dir,"Pertussis-calibration-",x,".trace"), sep = "")))
                              ids <- rep(x, length(trace))
                              # rnk <- rep(rank$distance[which(rank$id==x)]/max(rank$distance), length(trace))
                              rnk <- rep(rank$distance[which(rank$id==x)], length(trace))
                              trace <- data.frame(Trace=trace,id=ids, Time=1974+c(0:(length(unlist(trace))-1)), Score=rnk )
                              return(trace)
                            }
                      },
                    rnk=rank,
                    dir=folder,
                    cl=cl
                    )

stopCluster(cl = cl)
traces <- as.data.frame(do.call("rbind",traces))


trace <- filter(read.csv(paste0(folder,"Pertussis-calibration-",rank[1,]$id,".trace"), sep = ""))

trace <- data.frame(Trace=trace, Time=1974+c(0:(length(unlist(trace))-1)))

hyst <- read.csv("./input/reference_data.csv", header = FALSE, sep = "") # [1:21]

hyst <- data.frame(Trace=t(as.data.frame(hyst)), Time=1974+c(0:(length(unlist(hyst))-1)))

#traces <- traces[order(traces$Score, decreasing = FALSE),]


# s <- summary(traces$Score)
# s1 <- summary(traces$Score[-which(traces$Score>s[5])])
# cl <- makeForkCluster(nnodes = par_processors)
# idx_q3 <-parSapply(X=seq(from=s1[1],by=(s1[6]-s1[1])/5e4,to=s1[6]),
#                    FUN=function(x,traces){traces$id[which.min(abs(traces$Score-x))]},
#                    traces=traces,
#                    cl=cl)
# idx_q3 <-unique(idx_q3)
# # s2 <- summary(traces$Score[-c(which(traces$Score>s[5]), which(traces$Score < s1[5]))])
# s2 <- summary(traces$Score[which(traces$Score>s[5])])
# idx_q4_q3 <-parSapply(X=seq(from=s2[1],by=(s2[5]-s2[1])/1e5,to=s2[5]),
#                      FUN=function(x,traces){traces$id[which.min(abs(traces$Score-x))]},
#                      traces=traces,
#                      cl=cl)
# idx_q4_q3 <-unique(idx_q4_q3)
# s3 <- summary(traces$Score[which(traces$Score>s2[5])])
# idx_q4_q4_q3 <-parSapply(X=seq(from=s3[1],by=(s3[5]-s3[1])/1e4,to=s3[5]),
#                       FUN=function(x,traces){traces$id[which.min(abs(traces$Score-x))]},
#                       traces=traces,
#                       cl=cl)
# idx_q4_q4_q3 <-unique(idx_q4_q4_q3)
# 
# stopCluster(cl = cl)
# idx <- c(idx_q3,idx_q4_q3,idx_q4_q4_q3)
# trcs <-traces[which(traces$id %in% idx),]
# trcs<-trcs[order(trcs$Score, decreasing = TRUE),]
traces <- traces[order(traces$Score, decreasing = TRUE),]

traces <- traces[which(log(traces$Score)<=20),]

plt <- ggplot()+
  geom_line(data=traces,aes(y=Trace,group=id,col=(Score-min(Score))/(max(Score)-min(Score)),x=Time)) +
  geom_line(data=trace,aes(y=Trace,x=Time),col="black", size = 1.2, linetype = 1) +
  geom_line(data=hyst,aes(y=Trace,x=Time),col="orangered2", size = 1.2, linetype = 1) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold"))+
  labs(x="Years", y="Number of cases")+
  scale_colour_gradientn("Error",
                         colours = c("#FFFF00","#33CC33","#330066"),
                         values= c(0,.07,1),
                         breaks=c(0,1),
                         labels=c("min","max"))+
  scale_x_continuous(limits=c(1974,1994),breaks = seq(1974,1994,5) ) +
  coord_cartesian(ylim=c(0,4e4))

ggsave(plt, file="~/Plots/model_calibration.pdf", device="pdf", dpi="retina", width=10, height=7.5, units="in")


