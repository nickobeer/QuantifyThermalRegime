

library(zoo)
library(tidyverse)
library(lubridate)

#==== OVERVIEW ====
#   Data files needed. These can be obtained from primary sources also. See below in script. ====
#   "sweApril1data.csv"
#   "Lewiston.Air.csv"
#    Water data is from direct web query

#   Files generated:
#   Dataframe written to "anatone.2023.data.frame.csv"
#   "airregime.pdf"
#   "waterregime.pdf"
#   "Results.allyears.pdf"

#==== FUNCTIONS  all start with "z" ====
source("AnaFunctions.R")


#=== Anatone day lengths

#==== daylength data and fit to 3 par =====
dl <- read.csv("data/anatone.daylengths.csv")
daylength <- cbind.data.frame(day=1:365,time=NA)
doy <- 0
for(i in 2:13){for(j in 1:31){
  if(!is.na(dl[j,i])){
    doy <- doy + 1; 
    daylength$day[doy] <- doy;
    hold <- as.numeric(str_split(dl[j,i],":")[[1]]); 
    daylength$time[doy]<- hold[1] + hold[2]/60
  }
}}


#==== Water Flow & Temp. data ====
 # Flow data are not used in the current model

enddate <- "2022-09-30"
string3.0 <- "https://waterdata.usgs.gov/nwis/dv?cb_00010=on&cb_00060=on&format=rdb&site_no=13334300&legacy=&referred_module=sw&period=&begin_date=1960-10-01&end_date="
string3 <- paste0(string3.0,enddate)
data0 <- read_delim(string3,comment="#")[-1,] 
data0 <- data0 %>% rename(ana.flow='174089_00060_00003')  %>%
  rename(max='174090_00010_00001') %>%
  rename(min = '180456_00010_00002')  
data0$ana.flow <- as.numeric(data0$ana.flow)
data0$max <- as.numeric(data0$max)
data0$min <- as.numeric(data0$min) 
data0 <- data0 %>% mutate(Date = as.Date(datetime)) %>%
  mutate(year = year(Date)) %>%
  mutate(day = yday(Date)) %>%
  mutate(ana.temp=(max+min)/2) %>%
  dplyr::select(Date,year,day,ana.flow,ana.temp)

data01 <- cbind.data.frame(data0,zDOY.WY(data0[,c(2,3)])[,3:4]) %>%
  dplyr::select(Date,WY,WYdoy,ana.flow,ana.temp) %>% rename(year=WY) %>% rename(doy=WYdoy)

data01 <- data01[data01$year <= 2022,]
data02 <- cbind.data.frame(data01,ana.temp.fill = zFiller(data01$ana.temp))
data02$ana.flow.fill = zFiller(data02$ana.flow) # potentially allows filling across the year boundary

table(data02$year,!is.na(data02$ana.temp))

# eliminate years with more than 35 missing daily values
data03 <- data02[data02$doy != 366,]   # clip all years to 365 days
data03 <- data03[data03$year != 1961 & 
                   data03$year != 1969 &
                   data03$year != 1970 & 
                   data03$year != 1975 & 
                   data03$year != 1983 & 
                   data03$year != 1984 & 
                   data03$year != 1985 & 
                   data03$year != 1986 & 
                   data03$year != 1987, ] 
                   # data03$year != 1987 & 
                   # data03$year != 1992, ]
data03$ana.temp.raw=data03$ana.temp  # retain the original data values
data03$ana.temp=data03$ana.temp.fill # in case there were fills on smaller missing sequences, use these data
data03$ana.flow=data03$ana.flow.fill
data03 <- data03 %>% dplyr::select(-c(ana.flow.fill,ana.temp.fill))

#==== Air data ====

# Lewiston Air 
# Request was made at web site behind https://www.ncdc.noaa.gov/cdo-web/
# site is Lewiston,ID
# data set is 'Daily Summaries'
# data through a complete water year (ends 9-30)

# ONLY AIR and PRECIP on version 1
air0 <- read_csv("./Lewiston.Air.csv")

air0 <- as_tibble(air0) %>% mutate(Date=as.Date(DATE,format="%m/%d/%Y")) %>% 
  mutate(TavgComp = (TMAX + TMIN)/2 ) %>% mutate(Precip=PRCP) %>%
  dplyr::select(-c(STATION,NAME,DATE,TAVG,TMAX,TMIN,SNOW,SNWD)) 
air0$TavgComp <- zFiller(air0$TavgComp)
air0 <- air0[!is.na(air0$TavgComp),]
air0 <- air0 %>% mutate(year = year(Date)) %>% mutate(jul=yday(Date))%>% mutate(Tair = (TavgComp-32)*5/9  ) 
# Move everything to  Water Year time. Clip for 365 day year and end WY 2022

air1 <- cbind.data.frame(air0,zDOY.WY(cbind.data.frame(air0$year,air0$jul))[,3:4])
air1 <- air1[air1$WY >= 1961 & air1$WY <= 2022,]
air1 <- air1[air1$WYdoy <= 365,]
air1$Precip[is.na(air1$Precip)] <- 0 

air2 <- air1  %>% rename(doy=WYdoy) %>% dplyr::select(-year) %>% rename(year=WY)

df <- left_join(data03,air2 %>% dplyr::select(-Date),by=c("year","doy")) 

#==== Air regime metrics  ====
metrics <- NULL
for(y in unique(air2$year)){
  I <- air2$year == y
  temps <- df$ana.temp[df$year==y]
  air <- air2$Tair[I]
  springairA <-  mean(air[183:212],na.rm=TRUE)
  springairMar <-  mean(air[152:182],na.rm=TRUE)
  springairMay <-  mean(air[213:244],na.rm=TRUE)
  gMeanWater <- mean(temps,na.rm=TRUE)
  gMedianWater <- median(temps,na.rm=TRUE)
  summerairJA <- mean(air[274:335],na.rm=TRUE)
  winterairJF <- mean(air[93:151],na.rm=TRUE)
  winterairDJ <- mean(air[62:120],na.rm=TRUE)
  gMeanAir <- mean(air,na.rm=TRUE)
  gMedianAir <- median(air,na.rm=TRUE)
  
  metrics <- rbind.data.frame(metrics,
                              cbind.data.frame(year=y,gMeanAir,gMedianAir,gMeanWater,gMedianWater,springairA,summerairJA,winterairDJ,winterairJF), make.row.names= F)
}


#==== Snow Water Equivalent ====
# Can load file 

# or BUILD from sources wit the swe.grabber.standalone.R


if(1){
  # SWEPEAK load in with peak value and peak day
  # metricz <- left_join(metrics,swepeak %>% dplyr::select(c(year,mean,day)),by="year") %>% rename(swe=mean) %>% mutate(peakday=day+92) %>% dplyr::select(-"day")
  metricz <- left_join(metrics,swepeak %>% dplyr::select(c(year,mean)),by="year") %>% rename(swe=mean) 

}


#==== WAS Fit 5 param Air regime ====
# airP5.init <- c(11.64,11.04,164.3,20,0)  # initial air regime parameters
# I <- air2$year >= 1988
# Airfit5 <- optim(par=airP5.init,fn=zSineAir5fit,df=air2[I,],method="BFGS",control=list(trace=3))
# # This is the air regime: Modeled air temperature
# AirModelpar <- Airfit5$par
# names(AirModelpar) <- c("M","N","P","J","K")

# update from JFaulkner
airP5.init <- c(log(2), 11.64,11.04,164.3,20,0)  # initial air regime parameters, ADDED log(sigma)

I <- air2$year >= 1988
Airfit5 <- optim(par=airP5.init,fn=zSineAir5fit,df=air2[I,],method="BFGS",control=list(trace=3), hessian=TRUE)
# This is the air regime: Modeled air temperature
AirModelpar1 <- Airfit5$par
names(AirModelpar1) <- c("log.sigma", "M","N","P","J","K")
AirModelpar <- Airfit5$par[-1]
names(AirModelpar) <- c("M","N","P","J","K")

## JF Example usage. CI for sigma included: 
Airfit5.cov = solve(Airfit5$hessian)  #returns the inverse
Airmodelpar1.SE = sqrt(diag(Airfit5.cov))  #SE's of estimated parameters
Airfit5.summary = data.frame(AirModelpar1, getFitSummary(Airfit5))
# For returning CI of sigma after fitting sigma on the log scale
# The easy way is to do the following:
Airfit5.sigma = exp(Airfit5.summary$ParmEst[1])
Airfit5.sigma.CI = exp(c(Airfit5.summary$CI.low[1], Airfit5.summary$CI.up[1]))


Mair <- Mregimeair <- cbind.data.frame("doy"=1:365,AirRegime=zSineAir5(AirModelpar,1:365))

#=== Assemble df ====
usedf <- left_join(df,metricz,by="year")
usedf <- left_join(usedf,Mair,by="doy")

# usedf is a long timeseries. Not all years have SWE for fitting 
print(head(usedf))
print(length(unique(usedf$year)))
print(length(unique(usedf$year[!is.na(usedf$swe)])))

ANA.DF <- usedf[!is.na(usedf$swe),]
print(t(head(ANA.DF)))

#==== Fit 6 param Water regime ====

AnaModel6 <- optim(par=c(14,12,140,9,200,300), fn=zHiatusFit, 
                  data=ANA.DF$ana.temp,
                  method="BFGS",
                  control=list(maxit=100000,trace=2), hessian=TRUE) 

AnaModel6par <- AnaModel6$par
names(AnaModel6par) <- c("AA","BB","CC","DD","EE","FF")
Regime <- zSinFuncHiatusParVec(1:365,AnaModel6par)

# diagnostics
## JF Example usage. CI for sigma included: 
AnaModel6.cov = solve(AnaModel6$hessian)  #returns the inverse
AnaModel6par1.SE = sqrt(diag(AnaModel6.cov))  #SE's of estimated parameters
AnaModel6.summary = data.frame(AnaModel6par, getFitSummary(AnaModel6))
# For returning CI of sigma after fitting sigma on the log scale
# The easy way is to do the following:
AnaModel6.sigma = exp(AnaModel6.summary$ParmEst[1])
AnaModel6.sigma.CI = exp(c(AnaModel6.summary$CI.low[1], AnaModel6.summary$CI.up[1]))



#==== Fit 5 param Water regime ====

AnaModel5 <- optim(par=c(14,12,159.3,9,200,300), fn=zHiatusFit, 
                  lower=c(10,6,159.1,0,125,280), upper=c(16,15,159.5,10,250,350),
                  data=ANA.DF$ana.temp,
                  method="L-BFGS-B",
                  control=list(maxit=100000,trace=2)) 

AnaModel5par <- AnaModel5$par
names(AnaModel5par) <- c("AA","BB","CC","DD","EE","FF")

Regime5 <- zSinFuncHiatusParVec(1:365,AnaModel5par)

# The final data frame. Is written to file below after other results
ANA.DF <- left_join(ANA.DF,cbind.data.frame("doy"=1:365,"Regime"=Regime),by=c("doy"))

#==== Fit annual Water regime  ====
Eachyear6water <- NULL
anawateryears <- unique(df$year)
for(y in anawateryears){
  I <- df$year==y
  x <- df[I,]
  fit = optim(par=c(12,10,160,4,170,311), fn=zHiatusFit, data=x$ana.temp,method="L-BFGS-B",
              lower=c(8,6,140,0,125,280), upper=c(16,15,175,10,250,350),
              control=list(maxit=100000,trace=0)) 
  x2x <- zSinFuncHiatusParVec(1:365,fit$par)
  Eachyear6water <- rbind.data.frame(Eachyear6water,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}
names(Eachyear6water)[1:7] <- c("year","A","B","C","D","E","F")

#==== Fit 5 param annual Water regime  ====
Eachyear5water <- NULL
anawateryears <- unique(df$year)
for(y in anawateryears){
  I <- df$year==y
  x <- df[I,]
  fit = optim(par=c(12,10,159.2,4,170,311), fn=zHiatusFit, data=x$ana.temp,method="L-BFGS-B",
              lower=c(8,6,159.3,0,125,280), upper=c(16,15,159.4,10,250,350),
              control=list(maxit=100000,trace=0)) 
  x2x <- zSinFuncHiatusParVec(1:365,fit$par)
  Eachyear5water <- rbind.data.frame(Eachyear5water,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}
names(Eachyear5water)[1:7] <- c("year","A","B","C","D","E","F")

zplotfit(Eachyear5water$A,Eachyear6water$A)
zplotfit(Eachyear5water$B,Eachyear6water$B)
zplotfit(Eachyear5water$C,Eachyear6water$C)
zplotfit(Eachyear5water$D,Eachyear6water$D)
zplotfit(Eachyear5water$E,Eachyear6water$E)
zplotfit(Eachyear5water$F,Eachyear6water$F)

#==== create then export the trends (timseries fits) of the parameters ====
x <- Eachyear6water$year 
y <- Eachyear6water$A ; Atrendfit <- lm(y ~ x) 
y <- Eachyear6water$B ; Btrendfit <- lm(y ~ x) 
y <- Eachyear6water$C ; Ctrendfit <- lm(y ~ x) 
y <- Eachyear6water$D ; Dtrendfit <- lm(y ~ x) 
y <- Eachyear6water$E ; Etrendfit <- lm(y ~ x) 
y <- Eachyear6water$F ; Ftrendfit <- lm(y ~ x) 
rm(x); rm(y)


#==== Fit annual 5 Par Water regime  ====
#==== F is fixed
Eachyear5water <- NULL
anawateryears <- unique(df$year)
for(y in anawateryears){
  I <- df$year==y
  x <- df[I,]
  fit = optim(par=c(12,10,160,4,170,334), fn=zHiatusFit, data=x$ana.temp,method="L-BFGS-B",
              lower=c(8,6,140,0,125,333.9), upper=c(16,15,175,10,250,334.1),
              control=list(maxit=100000,trace=0)) 
  x2x <- zSinFuncHiatusParVec(1:365,fit$par)
  Eachyear5water <- rbind.data.frame(Eachyear5water,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}

names(Eachyear5water)[1:7] <- c("year","A","B","C","D","E","F")

# create then export the trends (timseries fits) of the parameters
x <- Eachyear5water$year 
y <- Eachyear5water$A ; Atrendfit <- lm(y ~ x) 
y <- Eachyear5water$B ; Btrendfit <- lm(y ~ x) 
y <- Eachyear5water$C ; Ctrendfit <- lm(y ~ x) 
y <- Eachyear5water$D ; Dtrendfit <- lm(y ~ x) 
y <- Eachyear5water$E ; Etrendfit <- lm(y ~ x) 
y <- Eachyear5water$F ; Ftrendfit <- lm(y ~ x) 
rm(x); rm(y)

#==== Fit annual 4 Par Water regime  ====
#==== F and C are fixed =
Eachyear4water <- NULL
anawateryears <- unique(df$year)
for(y in anawateryears){
  I <- df$year==y
  x <- df[I,]
  fit = optim(par=c(12,10,159,4,170,334), fn=zHiatusFit, data=x$ana.temp,method="L-BFGS-B",
              lower=c(8,6,158.9,0,125,333.9), upper=c(16,15,159.1,10,250,334.1),
              control=list(maxit=100000,trace=0)) 
  x2x <- zSinFuncHiatusParVec(1:365,fit$par)
  Eachyear4water <- rbind.data.frame(Eachyear4water,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}

names(Eachyear4water)[1:7] <- c("year","A","B","C","D","E","F")

# create then export the trends (timseries fits) of the parameters
x <- Eachyear4water$year 
y <- Eachyear4water$A ; Atrendfit <- lm(y ~ x) 
y <- Eachyear4water$B ; Btrendfit <- lm(y ~ x) 
y <- Eachyear4water$C ; Ctrendfit <- lm(y ~ x) 
y <- Eachyear4water$D ; Dtrendfit <- lm(y ~ x) 
y <- Eachyear4water$E ; Etrendfit <- lm(y ~ x) 
y <- Eachyear4water$F ; Ftrendfit <- lm(y ~ x) 
rm(x); rm(y)



print(AirModelpar)
print(AnaModel6par)
globalfixedpars <- c(AnaModel6par['CC'],
                     AnaModel6par['FF'],
                     AirModelpar['P'],
                     AirModelpar['J'],
                     AirModelpar['K'],
                     AnaModel6par['EE'])


#==== Fit annual air regime ====
Eachyear5air <- NULL
airyears <- unique(df$year)
for(y in airyears){
  I <- df$year==y
  x <- df[I,]
  # Remember to chope off the log(sigma) in the front
  fit <- optim(par=airP5.init[-1],fn=zSineAir5fit.v0,df=x,method="BFGS",control=list(trace=3))
  
  x2x <- zSineAir5(fit$par,days=1:365)
  Eachyear5air <- rbind.data.frame(Eachyear5air,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}
rm(x);rm(x2x)

names(Eachyear5air)[1:6] <- c("year","M","N","P","J","K")

#==== Fit Air regimes by year with fixed P,J, and K ====
Eachyear5airMN <- NULL
airyears <- unique(df$year)
for(y in airyears){
  I <- df$year==y
  x <- df[I,]
  fit <- optim(par=c(1,1),fn=zAnaFormAirfitMN, mydf=x,
               fixedpars=globalfixedpars,
               method="BFGS",control=list(maxit=1000,trace=3))
  
  x2x <- zSineAir2(fit$par,days=1:365,fixedpars=globalfixedpars)
  Eachyear5airMN <- rbind.data.frame(Eachyear5airMN,cbind.data.frame(year=y,t(fit$par),maxtemp=max(x2x),maxtempday=x$doy[x2x == max(x2x)]) )
}
rm(x);rm(x2x)


AnaAirMN <- optim(par=c(1,1),fn=zAnaFormAirfitMN, mydf=ANA.DF,
                  fixedpars=globalfixedpars,
                  method="BFGS",control=list(maxit=1000,trace=3))

#==== Fit Air regime with covariates ====
# AnaAir <- optim(par=c(1,1,1,1),fn=zAnaFormAirfit, mydf=ANA.DF,
#                 fixedpars=globalfixedpars,
#                 method="BFGS",control=list(maxit=1000,trace=3))
AnaAir <- optim(par=c(log(2),1,1,1,1),fn=zAnaFormAirfit.wsigma, mydf=ANA.DF,
                fixedpars=globalfixedpars,
                method="BFGS",control=list(maxit=1000,trace=3),hessian=TRUE)

  plot(ANA.DF$doy,ANA.DF$Tair,pch=16,col=rgb(10,10,10,3,NULL,25))
  for(y in unique(ANA.DF$year)){
    lines(1:365,zAnaFormAir(AnaAir$par[-1],mydf=ANA.DF[ANA.DF$year==y,],fixedpars=globalfixedpars),col="orange")
  }
  AnaAirpar <- AnaAir$par; 
  names(AnaAirpar) <- c("log.sigma","m0","m1","n0","n1")
  
  # Diagnostics
  ## JF Example usage. CI for sigma included: 
  AnaAir.cov = solve(AnaAir$hessian)  #returns the inverse
  AnaAirpar1.SE = sqrt(diag(AnaAir.cov))  #SE's of estimated parameters
  AnaAir.summary = data.frame(AnaAirpar, getFitSummary(AnaAir))
  # For returning CI of sigma after fitting sigma on the log scale
  # The easy way is to do the following:
  AnaAir.sigma = exp(AnaAir.summary$ParmEst[1])
  AnaAir.sigma.CI = exp(c(AnaAir.summary$CI.low[1], AnaAir.summary$CI.up[1]))
  
  
# Add the air parameters to the fixed set.
globalfixedpars <- c(globalfixedpars,AnaAirpar[-1])


alignside <- "right"

#==== Fit Daily Water model with 8 parameters  AKA type 1 ==== 
#                       1,  2 ,  3,  4,   5,    6,  7,   8 
waterpar.init    <-    c(9, .5,  9, .1, 170,  -.3,  .5,  .3)
AnaWater3type1 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=3,
                  method="BFGS",control=list(maxit=1000,trace=2))
AllPar3type1 <- AnaWater3type1$par
names(AllPar3type1) <- c("a0","a1","b0","b1","e0","e1","d","p0")

if(0){

  # AnaWater1 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=1,
  #                   method="BFGS",control=list(maxit=1000,trace=1))
  # AnaWater2 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=2,
  #                   method="BFGS",control=list(maxit=1000,trace=2))
  # 
  # AnaWater4 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=4,
  #                    method="BFGS",control=list(maxit=1000,trace=2))
  # AnaWater5 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=5,
  #                    method="BFGS",control=list(maxit=1000,trace=2))
  # AnaWater6 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=6,
  #                    method="BFGS",control=list(maxit=1000,trace=2))
  # AnaWater7 <- optim(par=waterpar.init[1:8], fn=zAnaFit, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=7,
  #                    method="BFGS",control=list(maxit=1000,trace=2))
  # # Three is preferred?
  # AllPar1 <- AnaWater1$par
  # AllPar2 <- AnaWater2$par
  # AllPar4 <- AnaWater4$par
  # AllPar5 <- AnaWater5$par
  # AllPar6 <- AnaWater6$par
  # names(AllPar1) <- c("a0","a1","b0","b1","e0","e1","d","p0")
  # names(AllPar2) <- c("a0","a1","b0","b1","e0","e1","d","p0")
  # names(AllPar4) <- c("a0","a1","b0","b1","e0","e1","d","p0")
  # names(AllPar5) <- c("a0","a1","b0","b1","e0","e1","d","p0")
  # names(AllPar6) <- c("a0","a1","b0","b1","e0","e1","d","p0")
  
}  # Block holding code for other lag values

# CORE RUN!
#==== Fit Daily Water model with 6 parameters AKA type=2 ==== 
# AnaFit2 does not try to fit E so the globalfixedpars value is used.
# Need to cherry pick the parameters sent to the fittter
# The "3" refers to the lag chosen

AnaWater3 <- optim(par=c(log(2),waterpar.init[c(1:4,7,8)]), fn=zAnaFit2.wsigma, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=3,
                   method="BFGS",control=list(maxit=1000,trace=2),hessian=TRUE)

AllPar3 <- AnaWater3$par[-1]
names(AllPar3) <- c("a0","a1","b0","b1","d","p0")
AllPar3.1 <- AnaWater3$par
names(AllPar3.1) <- c("log.sigma","a0","a1","b0","b1","d","p0")

AnaWater <- AnaWater3 ;useairlag <- 3 ; AllPar <- AllPar3
print(AllPar)

# Diagnostics
## JF Example usage. CI for sigma included: 
AnaWater.cov = solve(AnaWater$hessian)  #returns the inverse
AnaWaterpar1.SE = sqrt(diag(AnaWater.cov))  #SE's of estimated parameters
AnaWater.summary = data.frame(AllPar3.1, getFitSummary(AnaWater))
# For returning CI of sigma after fitting sigma on the log scale
# The easy way is to do the following:
AnaWater.sigma = exp(AnaWater.summary$ParmEst[1])
AnaWater.sigma.CI = exp(c(AnaWater.summary$CI.low[1], AnaWater.summary$CI.up[1]))



# For printing to manuscript table of coeffiients etc.:
print(Airfit5.summary)
print(AnaModel6.summary)
print(AnaWater.summary)
print(AnaAir.summary)




#==== Fit: Add term for autregressive correction. this is an experiment! ====
# # Adding another parameter for the autoregressive one day correction
# 
waterpar.init.a1    <- c(9, .5,  9, .1, 170,  -.3,  1,  .3, 0.3)

AnaWater3.a1 <- optim(par=waterpar.init.a1[c(1:4,7,8,9)], fn=zAnaFit2, daf=ANA.DF, fixedpars=globalfixedpars,useairlag=3,addAuto=TRUE,
                   method="BFGS",control=list(maxit=1000,trace=2))

AllPar3.a1 <- AnaWater3.a1$par
names(AllPar3.a1) <- c("a0","a1","b0","b1","d","p0","a1")














#==== Leave-one-out evaluation ====
# useairlag = 3 is preferred ?

if(1){  ##  slow. Do it only when the time for the 34 calibrations is no matter
  # set up for eiher the 8 parameter or 6 parameter version
  LOYOpars <- NULL
  years <- unique(ANA.DF$year)
  for(y in years){
    AnaModeloneout <- optim(par=waterpar.init[c(1:4,7,8)],fn=zAnaFit2,daf=ANA.DF[ANA.DF$year != y,],fixedpars=globalfixedpars,useairlag=3,
                            method="BFGS",control=list(maxit=10000,trace=2),addAuto=FALSE)
    LOYOpars <- rbind.data.frame(LOYOpars,c(y,AnaModeloneout$par))
  } ; cat("\n");
  # names(LOYOpars ) <- c("year","a0","a1","b0","b1","e0","e1","d","p0")
  names(LOYOpars ) <- c("year","a0","a1","b0","b1","d","p0")
} # end leaveoneout run

#==== leave one out analysis ====
if(1){ # leaveoneout analysis
  leaveoutresults <- NULL
  mod.obs.diffs.leaveone <- NULL
  mod.one.out.vals <- NULL
  useairlag <- 3

  for(y in years){
    legwords <- leglines <- legcols <- NULL
    AnaModeloneout <- LOYOpars[LOYOpars$year == y,][-1]
    
    K <- ANA.DF$year==y 
    # 8 param w/ zAnaFit
    # MD <- zAnaFit(AnaModeloneout,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag)
    # 6 param w/ zAnaFit2
    MD <- zAnaFit2(AnaModeloneout,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag)
    # Regime <- zAnaForm2(AnaModeloneout,ANA.DF[K,],fixedpars=globalfixedpars)
    airOb <- ANA.DF$Tair[K]
    Ob <- ANA.DF$ana.temp[K]
    x <- ANA.DF$doy[K]
    flow <- ANA.DF$ana.flow[K]
    airRegime <- zAnaFormAir(globalfixedpars[c("m0","m1","n0","n1")],ANA.DF[K,],fixedpars=globalfixedpars)
    Q <- ANA.DF$ana.flow[K]
    M1 <- zSinFuncHiatusParVec(1:365,Eachyear6water[Eachyear6water$year==y,2:7])
    air5fit <- zSineAir5(as.numeric(Eachyear5air[Eachyear5air$year==y,2:6]),days=1:365) 
    airR <- airOb - air5fit  # the one day residual
    air5R <- rollmean(airOb - air5fit,k=useairlag,align=alignside,fill=NA) # the smoothed residual using the useairlag
    Q <- ANA.DF$ana.flow[K]
    
    mad <- mean(abs(Ob-MD),na.rm=T)
    mre <- mean(Ob-MD,na.rm=TRUE)
    rmse <- sqrt(mean((MD-Ob)^2))
    diffs <- Ob - MD 
    cumerror <- sum(Ob - MD)
    cumsumerror <- cumsum(Ob-MD)
    worstcumsum <- max(abs(cumsumerror))
    R2 <- summary(lm(Ob ~ MD))$r.squared
    # meandiff7 <- mean(zrunavg(Ob-MD,7))
    meandiff7 <- mean(rollmean(Ob-MD,7,align=alignside,fill=NA),na.rm=T)
    
    # diff7days <- quantile((zrunavg(MD - Ob,7)),probs=c(0.5,0.75,0.80,0.825,0.85,0.875,0.90,1))
    diff7days <- quantile((rollmean(MD - Ob,k=7,align=alignside,fill=NA)),probs=c(0.5,0.75,0.80,0.825,0.85,0.875,0.90,1),na.rm=TRUE)
    summerMAD <- mean(abs(diffs[274:335]))
    winterMAD <- mean(abs(diffs[193:151]))
    leaveoutresults <- rbind.data.frame(leaveoutresults,
                                        cbind.data.frame(leaveoutyear=y,mad=mad,mre=mre,rmse=rmse,cumerror=cumerror,varExp=R2,
                                                         worstcumsum=worstcumsum,error7max=diff7days[["100%"]],
                                                         error7.90=diff7days[["90%"]],error7.80=diff7days[["80%"]],
                                                         error7.50=diff7days[["50%"]],summerMAD,winterMAD,meandiff7))
    mod.obs.diffs.leaveone <- c(mod.obs.diffs.leaveone,diffs)
    mod.one.out.vals <- c(mod.one.out.vals,MD)
    
  } # end for(year) loop
} # end leaveoneout analysis
print(leaveoutresults)
for(i in 2:14) print(paste0(names(leaveoutresults[i])," min:",min(leaveoutresults[,i])," max:",max(leaveoutresults[,i])))
OB <- ANA.DF$ana.temp
probss <- c(0,0.1,0.125,0.15,0.2,0.25,0.5,0.75,0.80,0.825,0.85,0.875,0.90,1)
# diff7days <- quantile((zrunavg(OB-mod.one.out.vals,7)),probs=probss)
# mediandiff7 <- median(zrunavg(OB-mod.one.out.vals,7))
# meandiff7 <- mean(zrunavg(OB-mod.one.out.vals,7))
# diff10days <- quantile(abs(zrunavg(OB-mod.one.out.vals,10)),probs=probss)
diff7days <- quantile((rollmean(OB-mod.one.out.vals,7,align=alignside,fill=NA)),probs=probss,na.rm=TRUE)
mediandiff7 <- median(rollmean(OB-mod.one.out.vals,7,align=alignside,fill=NA),na.rm=TRUE)
meandiff7 <- mean(rollmean(OB-mod.one.out.vals,7,align=alignside,fill=NA),na.rm=TRUE)
diff10days <- quantile(abs(rollmean(OB-mod.one.out.vals,10,align=alignside,fill=NA)),probs=probss,na.rm=TRUE)

print(paste(
  "Mean7day = ", meandiff7,
  "Median7day = ", mediandiff7,
  "Max7day = ", diff7days[["100%"]],
  "MAD = ",  mean(abs(OB-mod.one.out.vals),na.rm=TRUE),"\n",
  "RMSE = ", sqrt(mean((OB-mod.one.out.vals)^2,na.rm=TRUE)),"\n",
  "MRE = ", mean(OB-mod.one.out.vals,na.rm=TRUE),"\n",
  "R2 = ", summary(lm(OB ~ mod.one.out.vals))$r.squared,"\n"
))

print(summary(leaveoutresults))
plot(mod.obs.diffs.leaveone,pch=16,col=rgb(10,10,10,8,NULL,25)) ; abline(h=0,col="white",lwd=5); abline(h=0,col="darkblue",lwd=2)

#==== Put the annual regimes into the DF
ANA.DF2 <- cbind.data.frame(ANA.DF,annualRegime=NA,annualFit=NA,dailyFittype2 = NA,dailyFit=NA,dailyFittype2.a1=NA )

for(y in unique(ANA.DF2$year)){
  # type 1 is old method wit regression for 'E' is depricated
    K <- ANA.DF2$year==y 
    Regime <- zAnaForm2(AllPar3,ANA.DF2[K,],fixedpars=globalfixedpars)
    ANA.DF2$annualRegime[K] <- Regime
    ANA.DF2$annualFit[K] <- MD
    ANA.DF2$dailyFittype2[K] <- zAnaFit2(AllPar3,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag,addAuto=FALSE)
    ANA.DF2$dailyFit[K] <- zAnaFit(AllPar3type1,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag,addAuto=FALSE)
    ANA.DF2$dailyFittype2.a1[K] <- zAnaFit2(AllPar3.a1,ANA.DF[K,],fixedpars=globalfixedpars,retype="pred",useairlag=useairlag,addAuto=TRUE)
}

#==== Additive method for getting an AR term
# Additive method is needed where the AR term is found AFTER the weather term 
# necessary for our model to be correctly additive with differnt levels of weather detail.
dailymodelbase <- zAnaFit2(AllPar3,ANA.DF,fixedpars=globalfixedpars,retype="pred",useairlag=useairlag,addAuto=FALSE)
dailymodelerror <- dailymodelbase- ANA.DF$ana.temp
it <- arima(dailymodelerror,order=c(1,0,0))

predictedANA <- ANA.DF2$dailyFittype2 + it$coef[1]*c(0,lag(dailymodelerror,1)[-1])


print(args(z5))

#==== Generate results and write out files. ====
metaresultstype1 <- z5(pdff="Results.type1.pdf",type=1)
# for E fixed:
metaresults <- z5(pdff="Results.type2.pdf",type=2)
metaresults.a1 <- z5(pdff="Results.type2.a1.pdf",type=2,ar=TRUE)

zf(1,1); z5.prev.(pdff=NULL,usepars=AllPar3,useairlag=3,type=2);

zf(1,1); z5(pdff="Supplement.pdf",type=2,ar=T);
z5(years=2017, type=1)
z5(years=2017,type=2)
z5(years=2017,type=2,ar=T)

# usepars=AllPar, 
# myfixedpars=globalfixedpars,
# useairlag=3,
# years=unique(ANA.DF$year),
# addAuto=FALSE,
# myobs=NULL,
# myair=NULL

zf(2,2)
par(mar=c(4,4,1,1))
zplotfitDF(ANA.DF2,"annualRegime","ana.temp")
zplotfitDF(ANA.DF2,"annualFit","ana.temp")
zplotfitDF(ANA.DF2,"dailyFittype2","ana.temp")
zplotfitDF(ANA.DF2,"dailyFittype2.a1","ana.temp")


write_csv(ANA.DF2,"./anatone.dataframe.30.Sept.2025.csv")

#==== References to key R objects
print(AnaModel6par)
print(AirModelpar)
print(AllPar3.a1)
print(AllPar3)
print(AllPar3type1)


zf(2,3)
zplotfitDF(df=Eachyear6water,"year","A")
zplotfitDF(df=Eachyear6water,"year","B")
zplotfitDF(df=Eachyear6water,"year","C")
zplotfitDF(df=Eachyear6water,"year","D")
zplotfitDF(df=Eachyear6water,"year","E",pch=16)
zplotfitDF(df=Eachyear6water,"year","F")
zplotfitDF(df=Eachyear6water,"E","F")

zf(2,3)
zplotfitDF(df=Eachyear5air,"year","M")
zplotfitDF(df=Eachyear5air,"year","N")
zplotfitDF(df=Eachyear5air,"year","P")
zplotfitDF(df=Eachyear5air,"year","J")
zplotfitDF(df=Eachyear5air,"year","K")

zf(2,1)
zplotfitDF(df=metrics,"year","winterairDJ")
zplotfitDF(df=metrics,"year","winterairJF")

junkkkk <- metricz[metricz$year >= 1987,c(1,10)]  ; junkkkk$swe <- junkkkk$swe*2.54
zplotfitDF(df=junkkkk,"year","swe")

#==== PLOT: Thermal Regime and Daily Mean Temperatures ====
if(0){
mod <- zSinFuncHiatusParVec(paramvals = AnaModel6par)
# usedf is 52 years. ANA.DF only encompasses the SNOTEL/SWE available years
daymeans <- usedf %>% group_by(doy) %>% summarise(avg=mean(ana.temp),median=median(ana.temp))
regimedayerrormeans <- mod - daymeans$avg
regimedayerrormedians <- mod - daymeans$median
regimedayRMSE <- sqrt(mean(regimedayerrormeans^2))
regimeQs <- NULL
for(doy in 1:365){
  regimeQs <- rbind(regimeQs,quantile(usedf$ana.temp[usedf$doy==doy],probs=seq(0,1,by=0.05)))
}

for(type in c("","pdf")){
  if(type =="pdf")  pdf(file="waterregime.pdf",width=12,height=8)

    par(mar=c(5,3,1,3),cex.axis=1.5)
    zf(1,1)
    plot(1:365,mod,type="n",lwd=3,ylim=c(-1,25),ylab=bquote("Water Temperature "~degree*"C"),xlab="",axes=FALSE,cex.lab=1.5)
    coll <- rgb(10,9,9,3,NULL,25)
    for(i in 1:10){ 
      polygon(c(1:365,365:1),c(regimeQs[,i],rev(regimeQs[,22-i])),col=coll,border=NA)
    }
    box();zxWYdates()
    
    axis(2,las=1)
    axis(2,at=0:25,labels=F,las=1)
    lines(1:365,daymeans$avg,col="brown",lwd=3)
    # lines(1:365,mod,col="white",lwd=5)
    lines(1:365,mod,col="darkblue",lwd=5)
    lines(1:365,zSine(AnaModel6par),col="darkblue",lwd=2,lty=2)
    lines(1:365,regimedayerrormeans,col="orange",lwd=2)
    abline(h=c(-1,0,1),col="grey70")
    axis(4,at=c(-1,0,1),las=1)
    mtext(bquote("Regime - day of year mean "~degree*"C"),4,adj=0.3,line=1.7,cex=1.3)
    legend("topleft",bty="n",inset=c(0.04,0),legend=c("Water Temperature Regime and data\nAnatone, WA guage (1962-2022)"),adj=0,cex=1.5)
    legend("topleft",bty="n",inset=c(0.10,0.17),
           legend=c("Thermal Regime with 6 parameters","Day of year mean ","Regime - day of year mean","Data Cloud 1962-2022"),
           lwd=c(3,3,3,10),col=c("darkblue","brown","orange",coll),
           cex=1.25)
    cexx <- 1.5
    legend("bottomright",inset=c(0.25,0.15),bty="n",legend=c(
      paste("A =",round(AnaModel6par[1],2)),
      paste("B =",round(AnaModel6par[2],2)),
      paste("C =",round(AnaModel6par[3],0)),
      paste("D =",round(AnaModel6par[4],2)),
      paste("E =",round(AnaModel6par[5],0)),
      paste("F =",round(AnaModel6par[6],0))
      ),cex=cexx)
    print(paste("RMSE",sqrt(mean(regimedayerrormeans^2))))
    # arrows(1,AnaModel6par[1],300,AnaModel6par[1],code=3,angle=20,length=0.1,lwd=2)
    segments(10,AnaModel6par[1],320,AnaModel6par[1],lwd=2)
    rect(70,12,110,14,col="white",border=NA)
    text(90,AnaModel6par[1],"Mean: A",cex=cexx)
    arrows(1,(AnaModel6par[1]-AnaModel6par[2]),1,(AnaModel6par[1]+AnaModel6par[2]),code=3,angle=20,length=0.1,lwd=2)
    text(1,1.5,expression(Range: A %+-% B),adj=0.1,cex=cexx)
    arrows(365-AnaModel6par[3],25,365,25,code=3,angle=20,length=0.1,lwd=2)
    rect(280,24.5,320,25.5,col="white",border=NA)
    text(300,25,"Phase C",cex=cexx)
    junk <- (AnaModel6par[6]+AnaModel6par[5])/2
    junk2 <- (zSine(AnaModel6par,junk)+mod[junk])/2
    arrows(junk,mod[junk],junk,zSine(AnaModel6par,junk),code=3,angle=20,length=0.14,lwd=5,col="white")
    arrows(junk,mod[junk],junk,zSine(AnaModel6par,junk),code=3,angle=20,length=0.1,lwd=2)
    rect(205,junk2-0.5,240,junk2+0.5,col="white",border=NA)
    text(210,junk2,"Snow Effect: D",cex=cexx)
    segments(AnaModel6par[5],2.5,AnaModel6par[5],mod[AnaModel6par[5]],col="white",lwd=4)
    segments(AnaModel6par[5],2.5,AnaModel6par[5],mod[AnaModel6par[5]],lwd=2)
    segments(AnaModel6par[6],2.5,AnaModel6par[6],mod[AnaModel6par[6]],col="white",lwd=4)
    segments(AnaModel6par[6],2.5,AnaModel6par[6],mod[AnaModel6par[6]],lwd=2)
    text(AnaModel6par[5],2,"Begin Melt: E",cex=cexx,adj=0.35)
    text(AnaModel6par[6],2,"End Melt: F",cex=cexx)
    if(type =="pdf")  dev.off()
} # end plot fig type

print(paste("RMSE:",regimedayRMSE))
}

#==== PLOT: Air temperature at Lewiston and 5 parameter Air temperature regime ====
if(0){
  airmod <- zSineAir5(AirModelpar,1:365)
  dayairmeans <- usedf %>% group_by(doy) %>% summarise(avg=mean(Tair),median=median(Tair))
  regimeairdayerrormeans <- airmod - dayairmeans$avg
  regimeairdayerrormedians <- airmod - dayairmeans$median
  regimedayRMSE <- sqrt(mean(regimeairdayerrormeans^2))
  regimeairQs <- NULL
  for(doy in 1:365){
    regimeairQs <- rbind(regimeairQs,quantile(usedf$Tair[usedf$doy==doy],probs=seq(0,1,by=0.05)))
  }
  for(type in c("","pdf")){
    if(type =="pdf")  pdf(file="airregime.pdf",width=13,height=9)
      par(mar=c(3,5,1,3),cex.lab=1.5,cex.axis=1.3)
      plot(1:365,1:365,ylim=c(-10,30),xlab="",ylab=bquote("Air temperature "~degree*"C"),type="n",axes=F);box()
      coll <- rgb(10,9,9,3,NULL,25)
      for(i in 1:10){
        polygon(c(1:365,365:1),c(regimeairQs[,i],rev(regimeairQs[,22-i])),col=coll,border=NA)
      }
      zxWYdates()
      
        maxday <- c(1:length(airmod))[airmod==max(airmod)]
        lines(1:365,dayairmeans$avg,col="white",lwd=7)
        lines(1:365,dayairmeans$avg,col="brown",lwd=5)
        lines(1:365,airmod,lwd=9,col="white")
        lines(1:365,airmod,lwd=7,col="darkblue")
        AF5p <- AirModelpar
        lines(1:365,zSineAir5(c(AF5p[1],  AF5p[2], AF5p[3], 0 , 0 )),col="darkblue",lwd=2,lty=2)
        par(las=1)
        axis(2,at=c(seq(-5,30,by=5)))
        axis(2,at=c(-8:30),labels=FALSE)
        offsett <- 7
        lines(1:365,regimeairdayerrormeans-offsett,col="orange",lwd=2)
        abline(h=c(-1,0,1)-offsett,col="grey70")
        par(las=1)
        axis(4,at=c(-1,0,1)-offsett,labels=c("-1","0","1"),cex.axis=1.2)
        par(las=0)
        mtext(bquote("Regime - day of year mean "~degree*"C"),4,adj=0.3,line=1.7,cex=1.3)
        legend("topleft",bty="n",inset=c(0.04,0),legend=c("Air Temperature Regime and data\nLewiston, ID (1962-2022)"),adj=0,cex=1.5)
        
        legend("topleft",bty="n",inset=c(0.1,0.12),
               legend=c("Air temperature regime with 5 parameters","Air regime without asymmetry terms","Day of year mean","Regime - day of year mean","Data Cloud 1962-2022"),
               lwd=c(7,2,5,3,10),col=c("darkblue","darkblue","brown","orange","grey80"),lty=c(1,2,1,1,1),cex=1.2)
        
        cexx <- 1.4
        print(AF5p)
        legend("bottomright",inset=c(0.25,0.2),bty="n",legend=c(
          paste("M =",round(AF5p[1],2)),
          paste("N =",round(AF5p[2],2)),
          paste("P =",round(AF5p[3],1)),
          paste("J =",round(AF5p[4],1)),
          paste("K =",round(AF5p[5],1))
        ),cex=cexx)
        
        
        segments(1,AF5p[1],365,AF5p[1],lwd=2)
        rect(70,12,110,14,col="white",border=NA)
        text(90,AF5p[1],"Mean: M",cex=cexx)
        arrows(1,(AF5p[1]-AF5p[2]),1,(AF5p[1]+AF5p[2]),code=3,angle=20,length=0.1,lwd=2)
        text(1,0,expression(Range: M %+-% N),adj=0.1,cex=cexx)
        arrows(365-AF5p[3],27.5,365,27.5,code=3,angle=20,length=0.1,lwd=2)
        rect(233,26.7,273,28.2,col="white",border=NA)
        text(253,27.5," Phase: P",cex=cexx)
      
      
      segments(365-AF5p[5],5,365-AF5p[5],zSineAir5(c(AF5p),365-AF5p[5]),lwd=2)
      segments(365,5,365,29,lwd=2)
      arrows(365-AF5p[5],10,365,10,code=3,angle=20,length=0.1,lwd=2)
      text(265,10.25,adj=0,paste("Asymmetric Offset\nof Phase: K"),cex=cexx)
      
      arrows(365-AF5p[5]-91.5,zSineAir5(c(AF5p),365-AF5p[5]-91.5),365-AF5p[5]-91.5-AF5p[4],zSineAir5(c(AF5p),365-AF5p[5]-91.5),code=3,angle=20,length=0.1,lwd=4,col="white")
      arrows(365-AF5p[5]-91.5,zSineAir5(c(AF5p),365-AF5p[5]-91.5),365-AF5p[5]-91.5-AF5p[4],zSineAir5(c(AF5p),365-AF5p[5]-91.5),code=3,angle=20,length=0.1,lwd=2)
      # rect(200,21.2,262,22.8,col="white",border=NA)
      text(205,22,paste("Maximum Asymmetry: J"),cex=cexx)
      if(type =="pdf")  dev.off()
      
  } # End figure types: "" and "pdf"
  print(paste("RMSE:",regimedayRMSE))
}


airmod <- zSineAir5(AirModelpar,1:365)
dayairmeans <- usedf %>% group_by(doy) %>% summarise(avg=mean(Tair),median=median(Tair))
regimeairdayerrormeans <- airmod - dayairmeans$avg
regimeairdayerrormedians <- airmod - dayairmeans$median
regimedayRMSE <- sqrt(mean(regimeairdayerrormeans^2))
regimeairQs <- NULL
for(doy in 1:365){
  regimeairQs <- rbind(regimeairQs,quantile(usedf$Tair[usedf$doy==doy],probs=seq(0,1,by=0.05)))
}
#==== PLOT: COMBINED air and water ====

for(type in c("png")){   #  ""){    # },"pdf")){    c("pdf")){ #   
  if(type =="pdf")  pdf(file="allregime.pdf",width=11,height=13)
  if(type =="png")  png(file="allregime.png",width=800,height=1000)
  {
  par(mar=c(3,5,1,5),cex.lab=1.2,cex.axis=1.1)
  plot(1:365,1:365,ylim=c(-25,32),xlab="",ylab="",type="n",axes=F);box()
  par(mgp=c(3,0.5,0))
  zxWYdates() ; par(mgp=c(3,1,0))
  par(cex.lab=1.5,cex.axis=1.3)
  x <- c(daylength$time[274:365],daylength$time[1:273])
  lines(1:365,(x-min(x))*1 + min(x)-35,type="l",col="grey40",lwd=8)
  lines(1:365,(x-min(x))*1 + min(x)-35,type="l",col="white",lwd=1)
  # All ref lines:
  zvabline(v=c(82,264,301,314),f=c(0,0.12,0.85,0.49),coll="grey40",ltyy=3)
  text(c(82,264,301,314)-7, -24.5,c("","  Insolation Peak","Air Peak","Water Peak"),srt=90,las=0)
  
  legend("bottomleft",inset=c(0.15,0.06),bty="n",lty=1,col= "grey40",lwd=8, legend=c("Solar\nInsolation"),cex=1.5)
  legend("bottomleft",inset=c(0.15,0.06),bty="n",lty=1,col= "white",lwd=1, legend=c("Solar\nInsolation"),cex=1.5)
   y <- c(8,9,10,11,12,13,14,15,16)
  yat <- (y-min(x))*1 + min(x)-35
  axis(2,at = yat,labels=y,col.axis="grey40")
  mtext("Daylength (hours)",side=2,line=2.5,adj=0,las=3,cex=1.8,col="grey40")

  palegreen <- rgb(10,20,10,5,NULL,25)
  palegreen2 <- rgb(10,20,10,10,NULL,25)
  paleblue <- rgb(10,10,25,5,NULL,25)
  paleblue2 <- rgb(10,10,25,10,NULL,25)
  
  points(ANA.DF$doy,ANA.DF$Tair,col=palegreen,pch=16,cex=0.7);
  medz2 <- rep(NA,365)
  for(i in 1:365)medz2[i] <- median(ANA.DF$Tair[ANA.DF$doy==i],na.rm=T)
  mtext(expression(paste("Air temperature ",degree,"C")),2,line=2,col="darkgreen",cex=2,las=3)
  par(las=1)
  axis(2,at=c(-10,-5,0,5,10,15,20,25,30),col.axis="darkgreen")
  axis(2,at=seq(-10,30,by=1),labels=NA)
 

  maxday <- c(1:length(airmod))[airmod==max(airmod)]
  lines(1:365,dayairmeans$avg,col="white",lwd=7)
  lines(1:365,dayairmeans$avg,col="brown",lwd=5)
  # lines(1:365,dayairmeans$median,col="white",lwd=7)
  # lines(1:365,dayairmeans$median,col="orange",lwd=5)
  lines(1:365,airmod,lwd=9,col="white")
  lines(1:365,airmod,lwd=7,col="darkgreen")
  AF5p <- AirModelpar
  lines(1:365,zSineAir5(c(AF5p[1],  AF5p[2], AF5p[3], 0 , 0 )),col="darkgreen",lwd=2,lty=1)
 
  legend("topleft",bty="n",inset=c(0.0,0),legend=c("Air Thermal Regime and data\nLewiston, ID (1962-2022)"),adj=0,cex=1.4,text.col="darkgreen")
  legend("topleft",bty="n",inset=c(0.04,0.08),
         legend=c("Multi-year Air regime","Air regime w/o asymmetry","Day of year mean","Data Cloud 1962-2022"),
         lwd=c(7,2,5,10),col=c("darkgreen","darkgreen","brown",palegreen2),lty=c(1,1,1,1),cex=1)
  
  cexx <- 1.1
  print(AF5p)
  legend("topleft",inset=c(0.12,0.2),bty="n",legend=c(
    paste("M =",round(AF5p[1],2)),
    paste("N =",round(AF5p[2],2)),
    paste("P =",round(AF5p[3],1)),
    paste("J =",round(AF5p[4],1)),
    paste("K =",round(AF5p[5],1))
  ),cex=cexx,text.col="darkgreen")
  
  
  segments(-2,AF5p[1],70,AF5p[1],lwd=2)
  rect(70,12,110,14,col="white",border=NA)
  text(90,AF5p[1],"Mean: M",cex=cexx,col="darkgreen")
  arrows(1,(AF5p[1]-AF5p[2]),1,(AF5p[1]+AF5p[2]),code=3,angle=20,length=0.1,lwd=2)
  text(1,0,expression(Range: M %+-% N),adj=0.1,cex=cexx,col="darkgreen")
  arrows(365-AF5p[3],27.5,365,27.5,code=3,angle=20,length=0.1,lwd=2)
  rect(233,26.7,273,28.2,col="white",border=NA)
  text(253,27.5," Phase: P",cex=cexx,col="darkgreen")
  
  
  segments(365-AF5p[5],7,365-AF5p[5],zSineAir5(c(AF5p),365-AF5p[5]),lwd=2)
  segments(365,7,365,29,lwd=2)
  arrows(365-AF5p[5],10,365,10,code=3,angle=20,length=0.1,lwd=2)
  rect(272,9,320,12,col="white",border=NA)
  text(295,10.25,adj=0,paste("Asymmetry (K)\noffset to P"),cex=cexx,col="darkgreen")
  
  arrows(365-AF5p[5]-91.5,zSineAir5(c(AF5p),365-AF5p[5]-91.5),365-AF5p[5]-91.5-AF5p[4],zSineAir5(c(AF5p),365-AF5p[5]-91.5),code=3,angle=20,length=0.1,lwd=4,col="white")
  arrows(365-AF5p[5]-91.5,zSineAir5(c(AF5p),365-AF5p[5]-91.5),365-AF5p[5]-91.5-AF5p[4],zSineAir5(c(AF5p),365-AF5p[5]-91.5),code=3,angle=20,length=0.1,lwd=2)
  # rect(200,21.2,262,22.8,col="white",border=NA)
  text(190,21.5,paste("Maximum Asymmetry: J"),cex=cexx,col="darkgreen")
print(paste("RMSE:",regimedayRMSE))

  }

  {
mod <- zSinFuncHiatusParVec(paramvals = AnaModel6par)
# usedf is 52 years. ANA.DF only encompasses the SNOTEL/SWE available years
daymeans <- usedf %>% group_by(doy) %>% summarise(avg=mean(ana.temp),median=median(ana.temp))
regimedayerrormeans <- mod - daymeans$avg
regimedayerrormedians <- mod - daymeans$median
regimedayRMSE <- sqrt(mean(regimedayerrormeans^2))
regimeQs <- NULL

{

axis(4,at=c(-20,-15,-10,-5, 0),labels=c(0,5,10,15,20),col.axis="blue")
axis(4,at=seq(-20,0,by=1),labels=NA)
mtext(expression(paste("Water temperature ",degree,"C")),4,line=3,adj=0.2,cex=2,col="blue",las=3)
points(ANA.DF$doy,ANA.DF$ana.temp-20,col=paleblue,pch=16,cex=0.7)

lines(1:365,daymeans$avg-20,col="brown",lwd=3)
  # lines(1:365,mod,col="white",lwd=5)
  lines(1:365,mod-20,col="darkblue",lwd=5)
  lines(1:365,zSine(AnaModel6par)-20,col="darkblue",lwd=2,lty=1)

  legend("topright",bty="n",inset=c(0.0,0.0),legend=c("Water Thermal Regime and data\nAnatone, WA gauge (1962-2022)"),adj=0,cex=1.4,text.col="darkblue")
  rect(270,-17,330,-12,col="white",border=NA)
  legend("bottomright",border="white", bty="n", inset=c(0.10,0.16),
         legend=c("Multi-year Water Thermal Regime","Regime w/o snowmelt effect","Day of year mean ","Data Cloud 1962-2022"),
         lwd=c(7,2,5,10),col=c("darkblue","darkblue","brown",paleblue2),
         cex=1.1)
  cexx <- 1.1
  legend("bottomleft",inset=c(0.28,0.25),bty="n",legend=c(
    paste("A =",round(AnaModel6par[1],2)),
    paste("B =",round(AnaModel6par[2],2)),
    paste("C =",round(AnaModel6par[3],0)),
    paste("D =",round(AnaModel6par[4],2)),
    paste("E =",round(AnaModel6par[5],0)),
    paste("F =",round(AnaModel6par[6],0))
  ),cex=cexx,text.col="darkblue")
  print(paste("RMSE",sqrt(mean(regimedayerrormeans^2))))
  # arrows(1,AnaModel6par[1],300,AnaModel6par[1],code=3,angle=20,length=0.1,lwd=2)
  segments(-2,AnaModel6par[1]-20,70,AnaModel6par[1]-20,lwd=2)
  rect(38,12.3-20,70,13.7-20,col="white",border=NA)
  text(55,AnaModel6par[1]-20,"Mean: A",cex=cexx,col="darkblue")
  arrows(368,(AnaModel6par[1]-AnaModel6par[2])-20,368,(AnaModel6par[1]+AnaModel6par[2])-20,code=3,angle=20,length=0.1,lwd=2)
  text(370,0.5-20,expression(Range: A %+-% B),adj=1,cex=cexx,col="darkblue")
  arrows(365-AnaModel6par[3],25-20,365,25-20,code=3,angle=20,length=0.1,lwd=2)
  rect(280,24.5-20,320,25.5-20,col="white",border=NA)
  text(300,25-20,"Phase: C",cex=cexx,col="darkblue")
  junk <- (AnaModel6par[6]+AnaModel6par[5])/2
  junk2 <- (zSine(AnaModel6par,junk)+mod[junk])/2
  arrows(junk,mod[junk]-20,junk,zSine(AnaModel6par,junk)-20,code=3,angle=20,length=0.14,lwd=5,col="white")
  arrows(junk,mod[junk]-20,junk,zSine(AnaModel6par,junk)-20,code=3,angle=20,length=0.1,lwd=2)
  rect(205,junk2-0.5-20,240,junk2+0.5-20,col="white",border=NA)
  text(210,junk2-20,"Snow Effect: D",cex=cexx,col="darkblue")
  segments(AnaModel6par[5],2.5-20,AnaModel6par[5],mod[AnaModel6par[5]]-20,col="white",lwd=4)
  segments(AnaModel6par[5],2.5-20,AnaModel6par[5],mod[AnaModel6par[5]]-20,lwd=2)
  segments(AnaModel6par[6],16-20,AnaModel6par[6],mod[AnaModel6par[6]]-20,col="white",lwd=4)
  segments(AnaModel6par[6],16-20,AnaModel6par[6],mod[AnaModel6par[6]]-20,lwd=2)
  text(AnaModel6par[5]+10,1-20,"Begin Melt: E",cex=cexx,adj=0.35,col="darkblue")
  text(AnaModel6par[6],15.5-20,"End Melt: F",cex=cexx,col="darkblue",adj=0.5)
}
  

if(type =="pdf")  dev.off()
if(type =="png")  dev.off()
  } 
  }
# end plot fig type  pdf

print(paste("RMSE:",regimedayRMSE))


#==== Extra stuff. ====
# 
# junk <- zSineAir5(pv=c(13,12,270,20,12))
# begin <- 160
# end < 335
# mag <- 10
# dayz <- 1:365
# test <- dayz %% 366 >= begin & (dazy %% 366) <= end
# out <- rep(0,dayz)
# out2 <- out
# out2[test] <- out[test] - (mag/2 - mag/2 * sin((2 * pi)/(end - begin) * (day %% 365 - begin + (end - begin)/4)))[test]
# #  if(begin >= 0 & mag >= 0 & end > begin) {
# out[test] <- out2[test]
# 
# plot(1:365,junk,type="l")

for(i in 1:ncol(metricz) )print(paste(names(metricz)[i],signif(var(metricz[,i],na.rm=TRUE),3),signif(sqrt(var(metricz[,i],na.rm=TRUE)))))
sqrt(var(sweuseApril1$mean))
                                