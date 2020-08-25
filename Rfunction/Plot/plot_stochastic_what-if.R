library(ggplot2)
library(cowplot)
#######################
vaccination <-"80"
years<-"10"
folder_main <-paste0("~/results_stochastic_model_analysis_vaccination_",vaccination,"_years_",years,"/")
file_secondary <- "~/results_stochastic_model_analysis_vaccination_ref_years_10/Pertussis-analysys-1.trace"
#######################
main_tr <- read.csv(file=paste0(folder_main,"Pertussis-analysys-1.trace"), sep="")

secondary_tr <- read.csv(file=file_secondary, sep="")
#secondary_tr <- secondary_tr[1:(length(main_tr$Time)),]

n_run_m <- length(main_tr$Time)/length(unique(main_tr$Time))
n_run_s <- length(secondary_tr$Time)/length(unique(secondary_tr$Time))

col_names <- "(InfectCount_a){1}[0-9]{1}"
col_idxs <- c( which(names(main_tr) %in% grep(col_names, names(main_tr), value=T)) )
main_tr <- cbind(main_tr,data.frame(id=unlist(lapply(1:n_run_m, rep,times=44))))
secondary_tr <- cbind(secondary_tr,data.frame(id=unlist(lapply(1:n_run_s, rep,times=44))))
main_infects <- lapply(1:n_run_m,
                       function(x){
                         t <-main_tr[main_tr$id==x,]
                         tmp_inf <- rowSums(t[,col_idxs])
                         tmp_inf<-tmp_inf[-1]
                         diff <-c(tmp_inf[1], diff(tmp_inf,differences=1))
                         t<-t[-1,]
                         dt<-data.frame(Time=t$Time, Infects=diff, id=t$id)
                       })
secondary_infects <- lapply(1:n_run_s,
                            function(x){
                              t <-secondary_tr[secondary_tr$id==x,]
                              tmp_inf <- rowSums(t[,col_idxs])
                              tmp_inf<-tmp_inf[-1]
                              diff <-c(tmp_inf[1], diff(tmp_inf,differences=1))
                              t<-t[-1,]
                              dt<-data.frame(Time=t$Time, Infects=diff, id=t$id)
                            })
main_inf_trs<-do.call("rbind",main_infects)
secondary_inf_trs<-do.call("rbind",secondary_infects)

main_inf_trs$Time <- 1973+main_inf_trs$Time/365
secondary_inf_trs$Time <- 1973+secondary_inf_trs$Time/365

main_inf_trs$type <- "main_d"
secondary_inf_trs$type <- "a_secondary_d"


superdataframe<-rbind(secondary_inf_trs,main_inf_trs)

source("~/Rfunction/Plot/GeomSplitViolin.R")

plt<-ggplot(data = superdataframe[which(superdataframe$Time>2000),], aes(x=factor(Time),y=Infects,colour=type,fill=type))+
  geom_split_violin()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20,face="bold")) +
  labs( x="Years", y="Infects")+
  scale_fill_manual("Vaccination \npolicy",
                    values=alpha(c("main_d"="#000099","a_secondary_d"="#00CCFF"),.3),
                    labels=c("Real coverage",paste0("Fixed ", vaccination,"%")))+
  scale_colour_manual("Vaccination \npolicy",
                      values=c("main_d"="#000099","a_secondary_d"="#00CCFF"),
                      labels=c("Real coverage",paste0("Fixed ", vaccination,"%")))

ggsave(plot = plt,
       filename = paste0("~/Plots/ViolinPlot_p",vaccination,".pdf"),
       dpi = 400, width = 14, height = 8,device = "pdf")  


compare_ecdf <- ggplot(superdataframe[superdataframe$Time==max(superdataframe$Time),])+
  stat_ecdf(aes(Infects,group=type,col=type),geom = "step",size=1.1)+
  labs(y = "ECDF(Infects)", x="Infects")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20,face="bold"))+
  scale_colour_manual("Vaccination \npolicy",
                      values=c("main_d"="#000099","a_secondary_d"="#00CCFF"),
                      labels=c("Real coverage",paste0("Fixed ", vaccination,"%")))+
  xlim(1000,1500)
compare_ecdf

ggsave(plot = compare_ecdf,
       filename = paste0("~/Plots/ECDF_p",vaccination,".pdf"),
       dpi = 400,
       width = 14,
       height = 8,
       device = "pdf")

plt_col <-plot_grid(plt, compare_ecdf,
                    labels = c('a)', 'b)'),
                    ncol = 1
)

ggsave(plot = plt_col,
       filename = paste0("~/Plots/pv_", vaccination, "_", years, "years.pdf"),
       dpi = 400,
       width = 14,
       height = 16,
       device = "pdf")

