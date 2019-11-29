## Script authored by Michel Stefanello and modified by Santiago Botia 

### Functions.. 
fai_la_rotazione<-function(datiorari,na.rm=FALSE,zero_rot = TRUE,first_rot=TRUE,second_rot=FALSE){
  # datiorari: time,u,v,w,T
  ndati=length(datiorari[,1])
  uvw_old=datiorari[,c(1:3)]
  uvw_validi=subset(uvw_old,(is.na(uvw_old[,1])==FALSE | is.na(uvw_old[,2])==FALSE | is.na(uvw_old[,3])==FALSE),select=c(1,2,3))
  ndativalidi=length(uvw_validi[,1])
  u_old=uvw_old[,1]
  v_old=uvw_old[,2]
  w_old=uvw_old[,3]
  #   if (ndativalidi!=0){
  if (ndativalidi >1 ){
    # calcolo delle medie orarie
    umed=mean(u_old,na.rm=TRUE)
    vmed=mean(v_old,na.rm=TRUE)
    wmed=mean(w_old,na.rm=TRUE)
    # calcolo matrice di covarianza
    cov0R<-matrix(nrow=3,ncol=3)
    cov0R[1,1]=cov(u_old,u_old,use="na.or.complete")
    cov0R[1,2]=cov(u_old,v_old,use="na.or.complete")
    cov0R[1,3]=cov(u_old,w_old,use="na.or.complete")
    cov0R[2,1]=cov0R[1,2]
    cov0R[2,2]=cov(v_old,v_old,use="na.or.complete")
    cov0R[2,3]=cov(v_old,w_old,use="na.or.complete")
    cov0R[3,1]=cov0R[1,3]
    cov0R[3,2]=cov0R[2,3]
    cov0R[3,3]=cov(w_old,w_old,use="na.or.complete")
    # variabili logiche di controllo sulle rotazioni
    #      zero_rot=TRUE    # per fare la prima rotazione (<v>=0)
    #      first_rot=TRUE   # se vuoi la seconda rotazione (<w>=0)
    #      second_rot=TRUE  # per la terza rotazione (<v'w'>=0)
    # determinazione della matrice di rotazione
    rota=genero_matrice_rotazione(umed,vmed,wmed,cov0R,zero_rot,first_rot,second_rot,na.rm=FALSE)
    # calcolo del vettore vento ruotato  
    #      urot=rota[1,1]*u_old+rota[1,2]*v_old+rota[1,3]*w_old
    #      vrot=rota[2,1]*u_old+rota[2,2]*v_old+rota[2,3]*w_old
    #      wrot=rota[3,1]*u_old+rota[3,2]*v_old+rota[3,3]*w_old
    urot=rota[1,1]*uvw_old[,1]+rota[1,2]*uvw_old[,2]+rota[1,3]*uvw_old[,3]
    vrot=rota[2,1]*uvw_old[,1]+rota[2,2]*uvw_old[,2]+rota[2,3]*uvw_old[,3]
    wrot=rota[3,1]*uvw_old[,1]+rota[3,2]*uvw_old[,2]+rota[3,3]*uvw_old[,3]      
    datiorariruot<-data.frame(urot,vrot,wrot)
  } else  datiorariruot<-matrix(NA,nrow=ndati,ncol=3)
  structure(datiorariruot)
  #return(datiorariruot)
}

genero_matrice_rotazione<-function(uu,vv,ww,cov0R,zero_rot,first_rot,second_rot,na.rm=FALSE){
  # inizializzazione matrice unita'
  aa<-matrix(0,nrow=3,ncol=3)
  for(i in 1:3){
    for(j in 1:3){
      if (i==j) aa[i,j]<-1
    }
  }
  bb<-matrix(nrow=3,ncol=3)
  cc<-matrix(nrow=3,ncol=3)
  dd<-matrix(nrow=3,ncol=3)
  rota<-aa
  if (zero_rot | first_rot | second_rot){
    voriz=sqrt(uu*uu+vv*vv)
    # prima rotazione
    if(zero_rot){
      alfa=0.
      if(uu!=0. & vv!=0.) alfa=atan2(vv,uu)
      sinal=sin(alfa)
      cosal=cos(alfa)
      bb[1,1]=cosal
      bb[1,2]=sinal
      bb[1,3]=0.
      bb[2,1]=-sinal
      bb[2,2]=cosal
      bb[2,3]=0.
      bb[3,1]=0.
      bb[3,2]=0.
      bb[3,3]=1.
      rota=bb%*%rota
    }
    # seconda rotazione
    if(first_rot){
      gamma=0.
      if(ww!=0. & voriz>0.) gamma=atan2(ww,voriz)
      singa=sin(gamma)
      cosga=cos(gamma)      
      cc[1,1]=cosga
      cc[1,2]=0.
      cc[1,3]=singa
      cc[2,1]=0
      cc[2,2]=1.
      cc[2,3]=0.
      cc[3,1]=-singa
      cc[3,2]=0.
      cc[3,3]=cosga
      rota=cc%*%rota
    }
    # terza rotazione
    if(second_rot){
      fi=0.
      cov2R=rota%*%cov0R%*%t(rota)
      vp2=cov2R[2,2]
      vpwp=cov2R[2,3]
      wp2=cov2R[3,3]
      if (vpwp!=0. | (vp2-wp2)!=0.) fi=0.5*atan2((2.0*vpwp),(vp2-wp2))
      cosfi=cos(fi)
      sinfi=sin(fi)
      dd[1,1]=1.
      dd[1,2]=0.
      dd[1,3]=0.
      dd[2,1]=0.
      dd[2,2]=cosfi
      dd[2,3]=sinfi
      dd[3,1]=0.
      dd[3,2]=-sinfi
      dd[3,3]=cosfi
      rota=dd%*%rota
    }
  }
  return(rota)
}


# Merging starts here
path         = "/work/mj0143/b301034/Scrapbook_Analysis/ATTO_Observational_Analysis/ATTO_My_DataAnalysis/data/MicroMet_Data/ContinousMeasurements/81m/81m_2014-3-4-5-7-8-9/"  
#path_out     = "./Data/ATTO_81m_2014_03-09_u_v_w.1min2.txt"
readfiles_3m = list.files(path)
numdiv = 600
numtmp = 0
myresults <- data.frame()

for(kk in 1:length(readfiles_3m)){
  
  giorno_3m          = read.csv(paste(path,readfiles_3m[kk],sep="/"),na.strings = -999.99,sep = ",",header = TRUE)
  giorno_3m$T        = giorno_3m$T+273.16
  giorno_3m$datetime = as.POSIXct(giorno_3m$datetime, format="%Y-%m-%d %H:%M:%S")
  
  for (i in 1:24){
     numtmp = numtmp +1
     dati3m = giorno_3m[((i-1)*36000+1):(i*36000),]
     
     if(all(is.na(dati3m$u)) &  all(is.na(dati3m$v)) &  all(is.na(dati3m$T))) 
       {val_ini   = rep(NA,36000/numdiv)
        stat_1min = data.frame(data = rep(dati3m$time[1],36000/(numdiv)),
                               angolo = rep(NA,36000/(numdiv)),
                               thetamed = rep(NA,36000/(numdiv)), 
                               sigma_theta = rep(NA,36000/(numdiv)),
                               mean_u = rep(NA,36000/(numdiv)),
                               mean_v = rep(NA,36000/(numdiv)),
                               mean_w = rep(NA,36000/(numdiv)), 
                               mean_T = rep(NA,36000/(numdiv)),
                               sigma_u = rep(NA,36000/(numdiv)),
                               sigma_v = rep(NA,36000/(numdiv)), 
                               sigma_w = rep(NA,36000/(numdiv)),
                               sigma_T = rep(NA,36000/(numdiv)), 
                               uw = rep(NA,36000/(numdiv)),
                               vw = rep(NA,36000/(numdiv)), 
                               uT = rep(NA,36000/(numdiv)),
                               vT = rep(NA,36000/(numdiv)),
                               wT = rep(NA,36000/(numdiv)))}
     
     else{
       numero = 1:dim(dati3m)[1]
       val_ini     = rep(NA,36000/(numdiv))
       stat_1min   = data.frame(data = rep(dati3m$datetime[1],36000/(numdiv)),
       angolo      = rep(NA,36000/(numdiv)),
       thetamed    = rep(NA,36000/(numdiv)), 
       sigma_theta = rep(NA,36000/(numdiv)),
       mean_u      = rep(NA,36000/(numdiv)),
       mean_v      = rep(NA,36000/(numdiv)),
       mean_w      = rep(NA,36000/(numdiv)), 
       mean_T      = rep(NA,36000/(numdiv)),
       sigma_u     = rep(NA,36000/(numdiv)),
       sigma_v     = rep(NA,36000/(numdiv)), 
       sigma_w     = rep(NA,36000/(numdiv)),
       sigma_T     = rep(NA,36000/(numdiv)), 
       uw          = rep(NA,36000/(numdiv)),
       vw          = rep(NA,36000/(numdiv)),
       uT          = rep(NA,36000/(numdiv)),
       vT          = rep(NA,36000/(numdiv)),
       wT          = rep(NA,36000/(numdiv)))
       dati3m[,c("u","v","w")] = fai_la_rotazione(dati3m[,c("u","v","w")],na.rm=TRUE,zero_rot=TRUE,first_rot=TRUE,second_rot=FALSE) 
       dati3m$T    = dati3m$T-273.16
        
       for (kkk in 1:length(val_ini)){
           temp = seq(1+600*(kkk-1),(600*kkk),1)
           stat_1min$data[kkk]    = strptime(mean(dati3m$datetime[temp[1:2]]),"%Y-%m-%d %H:%M:%S",tz="UTC")
           stat_1min$mean_u[kkk]  = mean(dati3m$u[temp] ,na.rm = TRUE)
           stat_1min$mean_v[kkk]  = mean(dati3m$v[temp] ,na.rm = TRUE)
           stat_1min$mean_w[kkk]  = mean(dati3m$w[temp] ,na.rm = TRUE)
           stat_1min$mean_T[kkk]  = mean(dati3m$T[temp],na.rm = TRUE)
           stat_1min$sigma_u[kkk] = sd(dati3m$u[temp],na.rm = TRUE)
           stat_1min$sigma_v[kkk] = sd(dati3m$v[temp],na.rm = TRUE)
           stat_1min$sigma_w[kkk] = sd(dati3m$w[temp],na.rm = TRUE)
           stat_1min$sigma_T[kkk] = sd(dati3m$T[temp],na.rm = TRUE)
           stat_1min$uw[kkk]      = cov(dati3m$u[temp],dati3m$w[temp],use="na.or.complete")
           stat_1min$vw[kkk]      = cov(dati3m$v[temp],dati3m$w[temp],use="na.or.complete")
           stat_1min$uT[kkk]      = cov(dati3m$u[temp],dati3m$T[temp],use="na.or.complete")
           stat_1min$vT[kkk]      = cov(dati3m$v[temp],dati3m$T[temp],use="na.or.complete")
           stat_1min$wT[kkk]      = cov(dati3m$w[temp],dati3m$T[temp],use="na.or.complete")
           }
       myresults <- rbind(myresults, stat_1min)
       #filewrite = paste(path_out,sep="")
       #write.table(file=filewrite,stat_1min,row.names=FALSE, col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)    
     }
   }
}
write.table(myresults,file="./Data/ATTO_81m_2014_03-09_u_v_w.1min_3.txt",row.names=FALSE, col.names=FALSE,quote=FALSE,sep="\t")
