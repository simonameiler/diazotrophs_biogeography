
dat <- read.csv('~/Dropbox/Working/DIAZOTROPHS/data/Tang_2019_Supplement_CrossCheck.csv')

d <- data.frame(lat=dat$LATITUDE,lon=dat$LONGITUDE,depth=dat$DEPTH..m.,
                r=dat$Richelia.nifH.Gene..x106.copies.m.3.,
                t=dat$Trichodesmium.nifH.Gene..x106.copies.m.3.,
                ub=dat$UCYN.B.nifH.Gene..x106.copies.m.3.,
                ua=dat$UCYN.A1.nifH.Gene..x106.copies.m.3.)

d$latlon <- paste(d$lat,d$lon,sep='')
ids      <- unique(d$latlon)
N        <- length(ids)

df <- data.frame(lat=numeric(),lon=numeric(),richelia=numeric(),
                 trichodesmium=numeric(),ucyn_a=numeric(),ucyn_b=numeric())

for(i in 1:N){
  dd <- d[d$latlon==ids[i],]
  lat <- dd$lat[1]
  lon <- dd$lon[1]
  
  rmean   <- ifelse(sum(!is.na(dd$r))>0,  mean(dd$r, na.rm=TRUE), NA)
  tmean   <- ifelse(sum(!is.na(dd$t))>0,  mean(dd$t, na.rm=TRUE), NA)
  uamean  <- ifelse(sum(!is.na(dd$ua))>0, mean(dd$ua,na.rm=TRUE), NA)
  ubmean  <- ifelse(sum(!is.na(dd$ub))>0, mean(dd$ub,na.rm=TRUE), NA)
  
  dmax <- max(dd$depth,na.rm=TRUE)
  
  dftmp <- data.frame(lat=lat,lon=lon,
                      richelia=rmean*dmax,
                      trichodesmium=tmean*dmax,
                      ucyn_a=uamean*dmax,
                      ucyn_b=ubmean*dmax)
  df <- rbind(df,dftmp)
}
  
write.csv(file='~/dropbox/working/diazotrophs/data/correction_integrated.csv',df,row.names=FALSE) 

  