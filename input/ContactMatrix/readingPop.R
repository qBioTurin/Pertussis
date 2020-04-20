library(readxl)
library(readr)
 #Dati_Istat_Italia <- read_excel("Dati_Istat_Italia.xlsx", 
#                                  +     sheet = "Popolazione", n_max = 3)

 pop1974 <- read_csv("pop1974.csv")
 Pop<-data.frame(Età=pop1974$Età,Value=pop1974$Value)
 
 
 a0<-paste(0,"anni")
 a1<-paste(1:18,"anni")

 #######
 c(0,1,seq(5,65,5)) -> s1
 c(0,seq(4,69,5)) -> s2
 
 Age_class_all<-c(paste0(s1,"-",s2), "70++" )
 
 lapply(1:(length(s1)+1), function(x){
   
   if(x==(length(s1)+1)){
     a<- a2<-c(paste(19:100,"anni"),"100 anni e più")  
   }else{
     a<-paste(s1[x]:s2[x],"anni")
   }
   sum(Pop$Value[which(Pop$Età%in%a)])
   
 } ) -> l
 
 unlist(l)->Ni
 names(Ni)<-Age_class_all
 
########################################
contact <- read_csv("contact.csv") 

contact <- as.matrix(contact[,-1] )

contact_0<-contact[1,]*.2

contact_tmp<-rbind(contact_0,contact)
contact_tmp <- cbind(c(contact_0[1], contact_0),contact_tmp)
contact_tmp->contact
row.names(contact) =  Age_class_all
colnames(contact) =  Age_class_all



a0<-1
a1<-2:5
a2<-6:16

ages<-list(a0,a1,a2)

l_ages<-length(ages)

Mall<-matrix(0,ncol=l_ages,nrow = l_ages)

for(i in 1:l_ages)
{
  for (j in 1: l_ages)
  {
   if(length(Ni[ages[[i]]])>1) D<-diag(Ni[ages[[i]]])
    else D<-Ni[ages[[i]]]
    
    Mall[i,j] <- sum(D %*% contact[ages[[i]],ages[[j]]] )/sum(Ni[ages[[i]]]) 
    
  }
}


saveRDS(object = Mall, file = "ContactMatrix.RDS")
