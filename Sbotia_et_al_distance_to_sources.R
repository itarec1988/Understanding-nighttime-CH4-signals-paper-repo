# Script written by Christoph Gerbig, 2019
# Calculating distance to origin of CH4 enhancements

setwd("/mnt/lustre01/work/mj0143/b301034/Scrapbook_Analysis/ATTO_Observational_Analysis/ATTO_My_DataAnalysis/Paper1_UnderstandingCH4_ATTO_code_data_plots")

# This script should take the "ATTO_GHG_Inst_flux_soil_30min.csv" and generate the ATTO_30min_INSTall_Flux_CH4_Gerbig_dist2.csv

dat          <- read.csv("./Data/ATTO_30min_INSTall_Flux_CH4.csv") #Paper1_dataExcerpts/ATTO_30min_INSTall_Flux_CH4.csv")
dat$GradFlag <- toupper(dat$GradFlag)
ti           <- strptime(dat[,"LTime"], "%Y-%m-%d %H:%M")
dti          <- as.numeric(c(diff(ti),30))#in minutes

#head(dat[dat[,"GradFlag"]==TRUE,])

cum_wsp      <- cumsum(dat[,"WSp_73m"]*dti)

flg<-c(dat[2:length(dti),"GradFlag"]==TRUE&dat[1:(length(dti)-1),"GradFlag"]==FALSE,FALSE)
is20<-substring(ti,12,16)=="20:00"
is20s<-which(is20)

# Integrating distance in Km, based on the Wind Speed at 73 m
dist_ch4<-NULL
ids<-NULL
for(i in 1:(length(is20s)-1)){
  if(any(flg[is20s[i]:is20s[i+1]])) {
    flgi<-which(flg[is20s[i]:is20s[i+1]])[1]
    disti<-sum((dat[,"WSp_73m"]*dti*60)[is20s[i]+1:flgi])/1000 #distance in km
    dist_ch4<-c(dist_ch4,disti)
    ids<-c(ids,is20s[i]+flgi)
  } else {
    next
  }
}

length(dist_ch4)
dist_source73<-rep(NA,nrow(dat))
dist_source73[ids]<-dist_ch4
dat<-cbind(dat,dist_source73)

# Integrating distance in Km, based on the Wind Speed at 81 m
dist_ch4<-NULL
ids<-NULL
for(i in 1:(length(is20s)-1)){
  if(any(flg[is20s[i]:is20s[i+1]])) {
    flgi<-which(flg[is20s[i]:is20s[i+1]])[1]
    disti<-sum((dat[,"Mean_Windsp.m.s."]*dti*60)[is20s[i]+1:flgi])/1000 #distance in km
    dist_ch4<-c(dist_ch4,disti)
    ids<-c(ids,is20s[i]+flgi)
  } else {
    next
  }
}
dist_source81<-rep(NA,nrow(dat))
dist_source81[ids]<-dist_ch4
dat<-cbind(dat,dist_source81)

write.csv(dat,"./Data/ATTO_30min_INSTall_Flux_CH4_Gerbig_dist2.csv",row.names=FALSE)
#write.csv(dat,"TEST.csv",row.names = FALSE)
