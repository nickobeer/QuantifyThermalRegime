#===== SWE grabber =====
# Downloads SNOTEL data files
# Parses for snow water equivalent SWE metrics pertinent to the Snake River basin using
# Grande Rhonde, Powder, Burnt, Imnaha, Clearwtaer abd Salmon basins
# computes various metrics that are cataloged (see code for data.frames)
# see SNOTEL file documentation/headers  for the meaning of the metrics
# writes out the peak swe date and value for each year
# Pick the days Everday from  March 1  to  May 15
#                             Day  60  to  Day 135  # 76 days to search        # Leap year is ignored
library(tidyverse)
library(lubridate)

dir.create("swe")

# Download SNOTEL reports for every day and every year during March- mid May back to 1987
# 1987 is first year with full set of sites.
# server-side errors are NOT UNCOMMON. restart the download sequence at blocked day as needed 
# This is slow.
# To update the library, change years to a single year or other subset of all years.
years <- 1987:2024

# Download raw files
holdtime <- Sys.time()

#####==== VERY SLOW TO GRAB ALL  YEARS AND ALL DAYS =====
# for( i in 121:151){
for( i in 121){
  if(i<=31 & i >= 1){ mo <- "January"; da <- i } # no Feb 29 
  if(i<=59 & i >= 32){ mo <- "February"; da <- i-31 } # no Feb 29 
  if(i<=90 & i >= 60){ mo <- "March" ;  da <- i-59 }
  if(i<=120 & i >= 91){ mo <- "April" ;  da <- i-90 }
  if(i<=151 & i >= 121){ mo <- "May"   ;  da <- i-120 }
  if(i<= 181 & i >= 152){mo <- "June" ;  da <- i-151 }
  if(i<= 212 & i >= 182){mo <- "July" ;  da <- i-181 }
  if(i<= 243 & i >= 213){mo <- "August" ;  da <- i-212 }
  # Build the strings
  # https://wcc.sc.egov.usda.gov/reports/UpdateReport.html?textReport=Columbia+River+Basin&textRptKey=17&textFormat=SNOTEL+Snow%2FPrecipitation+Update+Report&StateList=Select+a+State&RegionList=17&SpecialList=Select+a+Special+Report&MonthList=April&DayList=8&YearList=1970&FormatList=N0&OutputFormatList=csv&textMonth=April&textDay=8&CompYearList=select+a+year
  string0 <- "https://wcc.sc.egov.usda.gov/reports/UpdateReport.html?textReport=Columbia+River+Basin&textRptKey=17&textFormat=SNOTEL+"
  # long? string02 <- "Snow%2FPrecipitation+Update+Report&StateList=Select+a+State&RegionList=17&SpecialList=Select+a+Special+Report"
  string02 <- "Snow%2FPrecipitation+Update+Report&RegionList=17"
  # Keep and inspect all files
  for(year in years){
  # long?  string1 <- paste("&MonthList=",mo,"&DayList=",da,"&YearList=",year,"&FormatList=N0&OutputFormatList=csv&textMonth=",mo,"&textDay=",da,"&CompYearList=select+a+year",sep="")
    string1 <- paste("&MonthList=",mo,"&DayList=",da,"&YearList=",year,"&FormatList=N0&OutputFormatList=csv&textMonth=",mo,"&textDay=",da,sep="")
    download.file(paste(string0,string02,string1,sep=""),method="libcurl",destfile=paste("swe/swe",year,".",mo,".",da,".txt",sep=""))
  }
  gc()
}

holdtime2 <- Sys.time()
print(holdtime2 -holdtime)

#==== Process loop for SWE ====

# Run through each raw report (sweYYYY.MM.DD.txt) and grab desired output and write to temporary files.
i.MoDa.lu <- NULL  # We want this to have a running sequence of all the dates across Month bounds.

for( i in 40:211){
# for(i in 91){
  if(i<=31 & i >= 1){ mo <- "January"; da <- i } # no Feb 29 
  if(i<=59 & i >= 32){ mo <- "February"; da <- i-31 } # no Feb 29 
  if(i<=90 & i >= 60){ mo <- "March" ;  da <- i-59 }
  if(i<=120 & i >= 91){ mo <- "April" ;  da <- i-90 }
  if(i<=151 & i >= 121){ mo <- "May"   ;  da <- i-120 }
  if(i<= 181 & i >= 152){mo <- "June" ;  da <- i-151 }
  if(i<= 212 & i >= 182){mo <- "July" ;  da <- i-181 }
  if(i<= 243 & i >= 213){mo <- "August" ;  da <- i-212 }
  # Read each year and begin parsing for relevant data
  swes <- swesum <- sweN <-  NULL
  x <- NULL
  for(year in years){
    tryfile <- paste("swe/swe",year,".",mo,".",da,".txt",sep="")
    if(file.exists(tryfile)){  x <- read.table(tryfile,sep=",",skip=74,header=TRUE)
          I  <- (x$Basin_name == "GRAND RONDE, POWDER, BURNT, IMNAHA" ) |
            (x$Basin_name == "CLEARWATER AND SALMON"& x$Lat < 45.7) 
            # | 
            # (x$Basin_name == "HENRYS FORK, TETON, WILLOW, BLACKFOOT, PORTNEUF" ) | 
            # (x$Basin_name == "WEISER, PAYETTE, BOISE" ) | 
            # (x$Basin_name == "RAFT, GOOSE, SALMON FALLS, BRUNEAU") | 
            # (x$Basin_name == "OWYHEE MALHEUR" ) |  
            # (x$Basin_name == "SNAKE ABOVE PALISADES")
          II <- x$Station_name != "Granite Creek" & x$Station_name != "Pine Creek Pass" & 
            x$Station_name != "Sedgwick Peak"  & x$Station_name != "Wilson Creek" & x$Station_name != "Jacks Peak" & x$Station_name != "Black Bear"  & 
            x$Station_name != "Crab Creek"  & x$Station_name != "Madison Plateau" & x$Station_name != "Whiskey Creek"
          xx <- x[I & II,][!is.na(x[I & II,1]),]
          swes <- rbind(swes,xx[,c(3,4,5,6,7,9,13)])
    }
    else { x <- NULL; print(paste0("File missing: ", tryfile)); swes <- NULL}
    }
  if(!is.null(swes)){
      # convert swe (inches) into cm !
      swes[,6] <- swes[,6] * 2.54
      samplesnotel <- x  # in case you want to keep a site
      swes <- swes %>% na.omit()  # remove NA rows Current R studio is having trouble RStudio.Version()$long_version ==  "2023.06.0+421"
      swes$Wteq_amt[swes$Wteq_amt == -999.0] <- NA  # NA for  missing values
      swes$Wteq_amt[swes$Wteq_amt <= -99.0] <- NA
      swes$Prec_wytd_amt[swes$Prec_wytd_amt <= -99.0] <- NA
      swesum <- swes %>% group_by(YYYYMMDD) %>% summarise(sum=sum(Wteq_amt,na.rm=TRUE)) 
      swemean <- swes %>% group_by(YYYYMMDD) %>% summarise(mean=mean(Wteq_amt,na.rm=TRUE))
      swemedian <- swes %>% group_by(YYYYMMDD) %>% summarise(median=median(Wteq_amt,na.rm=TRUE))
      precipsum <- swes %>% group_by(YYYYMMDD) %>% summarise(sum=sum(Prec_wytd_amt,na.rm=TRUE))
      precipmean <- swes %>% group_by(YYYYMMDD) %>% summarise(mean=mean(Prec_wytd_amt,na.rm=TRUE))
      precipmedian <- swes %>% group_by(YYYYMMDD) %>% summarise(median=median(Prec_wytd_amt,na.rm=TRUE))
      sweN  <- swes %>% group_by(YYYYMMDD) %>% summarise(N=n())
      for( ii in 1:nrow(sweN)){
        junk <- swes[swes$YYYYMMDD== sweN$YYYYMMDD[ii],]
        sweN$N[ii] <- length(junk$Wteq_amt[!is.na(junk$Wteq_amt)])  # compute counts in each year
      }
      # Write out the index file with the date you asked for
      swesum$year <- NA
      for(k in 1:nrow(swesum))swesum$year[k] <- as.numeric(substr(as.character(swesum[k,1]),1,4))
      
      # pp
      # Repair missing values
      # Write out the index file with the date you asked for
      assign(paste("sweuse",mo,da,sep=""),cbind(swesum,sweN,swemean,swemedian)[,c(1,3,5,2,7,9)])
      assign(paste("precipuse",mo,da,sep=""),cbind(precipsum,precipmean,precipmedian)[,c(1,2,4,6)])
      i.MoDa.lu <- rbind.data.frame(i.MoDa.lu,cbind.data.frame(i=i,moda=paste(mo,da,sep=""),mo=mo,da=da)) 
      cat(i,mo,da,"\n")
  } else {
    # assign(paste("sweuse",mo,da,sep=""),cbind(NA,NA,NA,NA)[,c(1,3,5,2,7,9)])
    # assign(paste("precipuse",mo,da,sep=""),cbind(NA,NA,NA)[,c(1,2,4,6)])
    # i.MoDa.lu <- rbind.data.frame(i.MoDa.lu,cbind.data.frame(i=i,moda=paste(mo,da,sep=""),mo=mo,da=da)) 
    # cat(i,mo,da,"\n")
  }

  }  # END  the loop for all the dates of interest

# Use xx to generate the table of locations
# unique(xx$Station_name)
# write.csv(file = "snotel.sites.csv",xx[,c(4,7,1,2,8)])

# diagnose i.MoDa.lu
# print(i.MoDa.lu)
# for(i in 60:212){
# for( i in 91){
#   print(objects(pattern=paste("sweuse",i.MoDa.lu[ i.MoDa.lu[,1] == i,2],sep="")))
# }

#=====
# Now open files in turn, by year and find the peak day / value and other mtricz  going through the days in sequence
junk <- NULL
swepeak <- NULL
swemeanarray <- swemedianarray <- matrix(NA,nrow=243,ncol=length(years))
for(year in years){
  
  for(i in 50:180){
    if(i == 50){now <- 0 ;   nowday <- 50}
    junk <- eval(as.name(paste("sweuse",i.MoDa.lu[ i.MoDa.lu[,1] == i,2],sep="")))
    wantedrow <- match(year,junk$year,nomatch=0)
    # cat(i," ",wantedrow," ")
    if(wantedrow > 0){
      test <- junk[wantedrow,]
      # now <- test$mean[1]
      # if(year==1990)print(test)
      # index the arrays with day of year
      thisdoy <- yday(ymd((test$YYYYMMDD)))
      # if(thisdoy != i) print(paste0("Evaluate an index mismatch here: i= ",i," and doy = ",thisdoy," year",year))
      
      # cat(test$mean[1], " ")
      # mismatch of day and index during leap years
      # try using only the index and ignoring doy
      # swemeanarray[thisdoy,match(year,years)] <- test$mean[1]
      # swemedianarray[thisdoy,match(year,years)] <- test$median[1]
      swemeanarray[i,match(year,years)] <- test$mean[1]
      swemedianarray[i,match(year,years)] <- test$median[1]
      if(i== 91)sweApril1mean <-  test$mean[1]
      if(i== 91)sweApril1median <-  test$median[1]
      if(!is.na(now) & !is.na(test$mean[1])){
        if(test$mean[1] > now){
          now <- test$mean[1]
          nowday <- i
          # print(paste(year,i,now,nowday))
        }
      }
    }
  }
  
  j <- eval(as.name(paste("sweuse",i.MoDa.lu[i.MoDa.lu[,1]==nowday,2],sep="")))[match(year,years),]
  out <- cbind.data.frame(year=year,date=i.MoDa.lu$moda[i.MoDa.lu[,1]==nowday],
                          day=nowday, # For precise day-of-year in leap years with date after Feb. 28 : add one
                          YYYMMDD=j[1,1],sum=j[1,4],N=j[1,3],mean=j[1,5],median=j[1,6],
                          Apr1mean=sweApril1mean,Apr1median=sweApril1median)
  # print(out)
  
  swepeak <- rbind.data.frame(swepeak,out)
}

#==== Replace swepeak with the simpler version withonly April 1 ====
sweApril <- NULL
# swemeanarray <- swemedianarray <- matrix(NA,nrow=243,ncol=length(years))

  i <- 91
  junk <- eval(as.name(paste("sweuse",i.MoDa.lu[ i.MoDa.lu[,1] == i,2],sep="")))
  sweApril <- junk # Simple. only April 1
write.csv(sweApril,file="sweApril.csv",row.names=FALSE)


#=== find some quantile dates in the swe trajectories ====
# this works if the full array of SWE have been scoured for the max day and max date etc.
swemeanfractions <- swemeanarray
swemedianfractions <- swemedianarray

for(year in years){
  j <- match(year,years)
  val1 <- swepeak$mean[swepeak$year == year]
  val2 <- swepeak$mean[swepeak$year == year]
  vec1<- swemeanarray[,j] 
  vec2<- swemedianarray[,j] 
 swemeanfractions[,j] <- vec1/rep(val1,length(vec1))
 swemedianfractions[,j] <- vec2/val2
 }

#==== Begin illustrations ==== 
makepdf <- TRUE
if(makepdf){pdf("meltprofiles.pdf",height=12,width=8)}
x <- swemeanarray
par(mar=c(3.5,4,1,1))
plot(x[,1],type="l",ylim=c(0,80),xlim=c(58,190),axes="F",xlab="",ylab="SWE")
box();axis(2);zxCaldates()
abline(v=91,col=rgb(10,10,10,10,NULL,25),lwd=3)
coll <- rainbow(40)
for(i in 1:ncol(x)){
  lines(x[,i],col=coll[i]) ; 
  # if(swepeak$year[i] == 2008 | swepeak$year[i] == 2011) lines(x[,i],col=coll[i],lwd=3) 
  text(59,x[60,i],years[i],adj=1,cex=0.5,col=coll[i])
  points(swepeak$day[i],swepeak$mean[i],col=coll[i],pch=16)
}
legend("topright",bty="n",inset=c(0.04,0.04),cex=1.2,legend=c("Points at peak day","April1 with heavy line"))
if(makepdf)dev.off()

#==== next illustrations ==== 
makepdf <- FALSE
if(makepdf){pdf("annualmeltprofiles.pdf",height=12,width=8)}
xarray <- swemeanarray
par(mar=c(3.5,4,1,4))
par(mfrow=c(2,2))
for(i in 1:ncol(xarray)){
    plot(xarray[,1],ylim=c(0,80),xlim=c(58,190),axes="F",xlab="",ylab="SWE",type="n")
    box();axis(2);zxCaldates()
    abline(v=91,col=rgb(10,10,10,10,NULL,25),lwd=3)
    
    lines(xarray[,i],col="black") ; 
    # if(swepeak$year[i] == 2008 | swepeak$year[i] == 2011) lines(xarray[,i],col=coll[i],lwd=3) 
    text(59,xarray[60,i],years[i],adj=1,cex=0.5,col="black")
    points(swepeak$day[i],swepeak$mean[i],col="black",pch=16)
    legend("topright",bty="n",inset=c(0.04,0.04),cex=1.2,legend=c("Points at peak day","April1 with heavy line"))
    
    whichyear <- years[i]
    use <- Eachyear6water %>% filter(year == whichyear) %>% dplyr::select(year,D,E,`F`)
    # use <- left_join(use,swemean,join_by(year)) 
    legend("topleft",bty="n",inset=c(0.05,0.05),legend=c(whichyear),cex=1.5)
    area <- rep(NA,nrow(use))
    par(new=TRUE)
    plot(0,0,xlim=c(120,350),ylim=c(0,10),type="n",xlab="DOY", ylab="",axes=FALSE)   
    axis(4)
    for(i in 1:nrow(use)){
      E=use$E[i]
      FF=use$F[i]
      D=use$D[i]
      x <- seq(E,FF,by=1)
      y <- zs(D,E,FF,i=x)
      area[i] <- zInt2(E,FF,D,E,FF)
      lines(x,y,col="blue",lwd=2)
    }
}
if(makepdf)dev.off()

#==== Quantiles grab ====
swepeak$q90 <- swepeak$q80 <- swepeak$q70 <- swepeak$q60 <- swepeak$q50 <- 
  swepeak$q40 <- swepeak$q30 <- swepeak$q20 <- swepeak$q10 <- NA
x <- swemeanfractions
plot(x[,1],type="l",ylim=c(0,1),xlim=c(58,190),axes="F",xlab="",ylab="SWE")
box();axis(2);zxCaldates()
coll <- terrain.colors(40)
for(i in 1:ncol(x)){
  lines(x[,i],col=coll[i]) ; if(swepeak$year[i] == 2008 | swepeak$year[i] == 2011) lines(x[,i],col=coll[i],lwd=3) 
  text(59,x[60,i],years[i],adj=1,cex=0.7,col=coll[i])
  points(swepeak$day[i],swepeak$mean[i],col=coll[i],pch=16)
  for(j in 1:nrow(x)){
    if(!is.na(x[j,i])){
      if(x[j,i]>= 0.9 & j > swepeak$day[i])swepeak$q90[i] <- j + 92
      if(x[j,i]>= 0.8 & j > swepeak$day[i])swepeak$q80[i] <- j  + 92
      if(x[j,i]>= 0.7 & j > swepeak$day[i])swepeak$q70[i] <- j  + 92
      if(x[j,i]>= 0.6 & j > swepeak$day[i])swepeak$q60[i] <- j  + 92
      if(x[j,i]>= 0.5 & j > swepeak$day[i])swepeak$q50[i] <- j  + 92
      if(x[j,i]>= 0.4 & j > swepeak$day[i])swepeak$q40[i] <- j  + 92
      if(x[j,i]>= 0.3 & j > swepeak$day[i])swepeak$q30[i] <- j  + 92
      if(x[j,i]>= 0.2 & j > swepeak$day[i])swepeak$q20[i] <- j  + 92
      if(x[j,i]>= 0.1 & j > swepeak$day[i])swepeak$q10[i] <- j  + 92
    }
  }
}

plot((swepeak$q10-swepeak$q90),swepeak$mean)

x <- swemeanarray
plot(x[,1],type="l",ylim=c(0,80),xlim=c(58,190),axes="F",xlab="",ylab="SWE")
box();axis(2);zxCaldates()
coll <- terrain.colors(40)
for(i in 1:ncol(x)){
  lines(x[,i],col=coll[i]); if(swepeak$year[i] == 2008 | swepeak$year[i] == 2011) lines(x[,i],col=coll[i],lwd=3) 
  text(59,x[60,i],years[i],adj=1,cex=0.7,col=coll[i])
  #points(swepeak$day[i],swepeak$mean[i],col=coll[i],pch=16)
  points(swepeak$q90[i]-92,swepeak$mean[i]*0.9,col=coll[i],pch=16)
  points(swepeak$q10[i]-92,swepeak$mean[i]*0.1,col=coll[i],pch=16)
}

write.table(swepeak,file="swedetails.csv",row.names=FALSE,sep=",")
junk <- left_join(Eachyear6water,swepeak,by=join_by(year))
