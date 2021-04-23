
#dat <- read.csv('~/Dropbox/Working/DIAZOTROPHS/data/Tang_2019_Supplement_CrossCheck.csv')
dat <- read.csv('~/Dropbox/Working/DIAZOTROPHS/data/Tang_2019_CMAP_v4.csv')

# d <- data.frame(lat=dat$LATITUDE,lon=dat$LONGITUDE,depth=dat$DEPTH..m.,
#                 r=dat$Richelia.nifH.Gene..x106.copies.m.3.,
#                 t=dat$Trichodesmium.nifH.Gene..x106.copies.m.3.,
#                 ub=dat$UCYN.B.nifH.Gene..x106.copies.m.3.,
#                 ua=dat$UCYN.A1.nifH.Gene..x106.copies.m.3.,
#                 ref=dat$SOURCE..Data)

d <- data.frame(lat=dat$LATITUDE,lon=dat$LONGITUDE,depth=dat$DEPTH..m.,
                r=dat$Richelia.nifH.Gene..copies.L.1.,
                t=dat$Trichodesmium.nifH.Gene..copies.L.1.,
                ub=dat$UCYN.B.nifH.Gene..copies.L.1.,
                ua=as.numeric(dat$UCYN.A1.nifH.Gene..copies.L.1.),
                ref=dat$SOURCE..Data)

d$latlon <- paste(d$lat,d$lon,sep='')
ids      <- unique(d$latlon)
N        <- length(ids)

df <- data.frame(lat=numeric(),lon=numeric(),richelia=numeric(),
                 trichodesmium=numeric(),ucyn_a=numeric(),ucyn_b=numeric(),ref=character())

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
                      ucyn_b=ubmean*dmax,
                      ref=dd$ref[1])
  df <- rbind(df,dftmp)
}

#dff <- df[df$ref=="Goebel et al. (2010), doi:10.1111/j.1462-2920.2010.02303.x",]
#plot(dff$lat,dff$ucyn_a)

write.csv(file='~/dropbox/working/diazotrophs/data/correction_integrated_04_23_2021.csv',df,row.names=FALSE) 

  