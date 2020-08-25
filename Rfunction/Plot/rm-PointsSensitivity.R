
load("./results_sensitivity_analysis/prcc_Pertussis-sensitivity.RData")
load("./results_sensitivity_analysis/parms_prcc_Pertussis-sensitivity.RData")

# Plot PRCC prob1= boost dopo infezione, prob2=infettare susceptible, prob3= infettare recovered in l1
sensitivity.3 <- read.csv("./results_sensitivity_analysis/Pertussis-sensitivity-3.trace", sep="")
# Get the parameter names and the total number of parameters
names_param=colnames(sensitivity.3)[c(1,2,3,26,27,50,51,52,74,75,76)+1]
names_param=paste0("init_",names_param)
names_param=c("prob_boost","prob_infectionS","prob_infectionR_l1",names_param)
n_params = length(names_param)
time <- 0:21


filesAll <- list.files('.',pattern='*.trace')

traceAll <- lapply(filesAll, function(x) {
  #IDsim<- -readr::parse_number(x)
  IDsim<- as.numeric(gsub(pattern="(.trace)", gsub(pattern="([[:graph:]]+(-){1})", x=x, replacement=""), replacement="") )
  data.frame(read.csv(paste(x), sep="") ,ID=IDsim)
} )

dati<-ldply(traceAll, data.frame)
indexes.s <- c( which(colnames(dati) %in% grep("Is_*", colnames(dati), value=T)) )
indexes.p <- c( which(colnames(dati) %in% grep("Ip_*", colnames(dati), value=T)) )


DatiPRCC<-dati[,c("Time","InfectCount_a1","InfectCount_a2","InfectCount_a3","ID")]
ID= unique(DatiPRCC[,"ID"])

dati.Infect<-lapply(ID ,function(x) {
  dati_Infect<-DatiPRCC[which(DatiPRCC[,"ID"]==x),-c(1,5)]
  dati_Infect<-rowSums(dati_Infect)[-1]
  dati_Infect<-c(dati_Infect[1],diff(dati_Infect,differences = 1))
  data.frame(Time=DatiPRCC[which(DatiPRCC[,"ID"]==x),"Time"][-1],InfectsCount=dati_Infect,IDpar=x)
}
)
dati.Infect<-ldply(dati.Infect, data.frame)


DatiPRCC.subset<-data.frame(Time=dati.Infect$Time,parms[dati.Infect[,3], ],OutPut=dati.Infect[,2])


##### plot con gli errori

reference<-c(7413,10786,18354,8076,12582,18142,14170,7335,16599,27146,11648,14607,16770,31244,8764,6183,15256,18845,6701,4279,13735)

DatiPRCC<-dati[,c("Time","InfectCount_a1","InfectCount_a2","InfectCount_a3","ID")]
ID= unique(DatiPRCC[,"ID"])



dati.Infect.error<-lapply(ID ,function(x) {
  dati_Infect<-DatiPRCC[which(DatiPRCC[,"ID"]==x),-c(1,5)]
  dati_Infect<-rowSums(dati_Infect)[-1]
  dati_Infect<-c(dati_Infect[1],diff(dati_Infect,differences = 1))
  ret <- sum(( dati_Infect - reference )^2 )
  data.frame(InfectsCount=ret,IDpar=x )
  #data.frame(t=0:20,Infects=dati_Infect,InfectsCount=rep(ret,length(dati_Infect)),IDpar=rep(x,length(dati_Infect)) )
}
)

dati.Infect.errors<-ldply(dati.Infect.error, data.frame)
colnames(parms)<-names_param

DatiError.subset<-data.frame(parms[dati.Infect.errors[,2], c("prob_infectionS","prob_infectionR_l1","init_S_a2") ],Error=(dati.Infect.errors[,1] - min(dati.Infect.errors[,1]))/max(dati.Infect.errors[,1])  )

DatiError.subset<-DatiError.subset[order(DatiError.subset$Error,decreasing = T),]

DatiError.subset<-DatiError.subset[DatiError.subset$prob_infectionS<0.25,]
DatiError.subset$Error<- (DatiError.subset$Error - min(DatiError.subset$Error))/max(DatiError.subset$Error)

ggplot(DatiError.subset,aes(x=prob_infectionS,y=prob_infectionR_l1,col=Error))+
  geom_point(size=2.5)+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20,face="bold")) +
  labs(title="")+
  scale_colour_gradientn(colours = c("black","deepskyblue2","cyan"),
                         values= c(0,.000001,1),
                         breaks=c(0,.999),
                         labels=c("min","max"))