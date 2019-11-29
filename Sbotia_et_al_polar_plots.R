# SBotia et al., 2019
## Nighttime methane signals at ATTO
## Previuos requirements: openair, lubridate and ggplot2

library(openair)
library(lubridate)
library(ggplot2)
library(dplyr)

path      = './Data/ATTO_GHG_Inst_flux_soil_30min.csv' 
fname_out = './Figures/'
data30min2013_2018v3_soil = read.csv(path)
names(data30min2013_2018v3_soil)[1] <- 'LocalTime'
df.79m2 <- data.frame(data30min2013_2018v3_soil$LocalTime,data30min2013_2018v3_soil$Mean_Windsp.m.s.,data30min2013_2018v3_soil$Wind.Direc.corrg.deg.,data30min2013_2018v3_soil$grad)
names(df.79m2) <- c('date','ws','wd','grad79m')
df.79m2$date   = ymd_hms(df.79m2$date)

# Figure 7. Annulus
png(file = paste(fname_out,'ATTO_annulus_diurnalLT.png',sep=""),height= 1600, width = 2000,res=300)
a <- polarAnnulus(df.79m2, poll = "grad79m", period = "hour",k=8,auto.text = FALSE,key.header = "[ppb]",key.footer = "Mean CH4 \n Gradient",force.positive = FALSE,width="fat",limits=c(-5,10),fontsize=10,key.position = "bottom")
dev.off()

png(file = paste(fname_out,'ATTO_annulus_seasonLT.png',sep=""),height= 1600, width = 2000,res=300)
b <- polarAnnulus(df.79m2, poll = "grad79m", period = "season",auto.text = FALSE,key.header = "[ppb]",key.footer = "Mean CH4 \n Gradient",force.positive = FALSE,width="fat",limits=c(-5,10),fontsize=10,key.position = "bottom")
dev.off()

png(file = paste(fname_out,'ATTO_annulus_season-diurnal-LT.png',sep=""),height= 1600, width = 2000,res=300)
print(a, split = c(1, 1, 2, 1))
print(b, split = c(2, 1, 2, 1), newpage = FALSE)
dev.off()

# Figure 8. Polar CPF
names(data30min2013_2018v3_soil)[148] = "Mean Wind Sp. 81m"
png(file = paste(fname_out,'ATTO_polar_percentiles_88-100_noFooter2.png',sep=""),height= 1600, width = 2000,res=300)
polarPlot(data30min2013_2018v3_soil,x='Mean Wind Sp. 81m',wd='Wind.Direc.corrg.deg.',pollutant = 'grad',statistic = 'cpf',percentile=c(88.33,100),force.positive = TRUE,
          exclude.missing = TRUE,upper=10,key.header = "CPF Prob.",key.footer = " ",limits=c(0,1),cols = c("grey", "tomato", "forestgreen"),grid.line=2,par.settings=list(fontsize=list(text=18)))
dev.off()

# Appendix 1 and 2. Wind Rose diurnal and monthly
png(file = paste(fname_out,'ATTO_WDir_diurnalLT.png',sep=""),height= 1600, width = 2000,res=300)
windRose(df.79m2,type='hour',layout = c(6,4), annotate = c(" ", " "),max.freq = 20,grid.line = 10,breaks = c(0, 1, 4, 7,30),fontsize=8)
dev.off()

png(file = paste(fname_out,'ATTO_WDir_monthlyLT.png',sep=""),height= 1600, width = 2000,res=300)
windRose(df.79m2,type='month',layout = c(6,2), annotate = c(" ", " "),max.freq = 20,grid.line = 10,breaks = c(0, 1, 4, 7,30),fontsize=8)
dev.off()


