library(R.matlab)
library(fields)
library(maps)
library(mapdata)
library(viridis)
data(world2HiresMapEnv)
data(world2LoresMapEnv)

setwd('d:/dropbox/working/diazotrophs/tang-cassar-2019/')

D <- readMat('Diazotrophs_RF_global_monthly.mat')  #read the full dataset

r  <- D$Richelia.RF.global.monthly  #extract individual groups
t  <- D$Tricho.RF.global.monthly
ua <- D$UCYNA.RF.global.monthly
ub <- D$UCYNB.RF.global.monthly

r[r=='1'] <- 0    #set 1s to true zeros
t[t=='1'] <- 0
ua[ua=='1'] <- 0
ub[ub=='1'] <- 0

lon <- seq(-180,180,length.out=180)
lat <- seq(-90,90,length.out=90)

#########################################################
## NIFH ANALYSIS ########################################
#########################################################
rmean  <- apply(r,c(1,2),mean)
tmean  <- apply(t,c(1,2),mean)
uamean <- apply(ua,c(1,2),mean)
ubmean <- apply(ub,c(1,2),mean)

rsd  <- apply(r,c(1,2),sd)
tsd  <- apply(t,c(1,2),sd)
uasd <- apply(ua,c(1,2),sd)
ubsd <- apply(ub,c(1,2),sd)

rrng   <- apply(r,c(1,2),function(x) diff(range(x)))
trng   <- apply(t,c(1,2),function(x) diff(range(x)))
uarng  <- apply(ua,c(1,2),function(x) diff(range(x)))
ubrng  <- apply(ub,c(1,2),function(x) diff(range(x)))

##--PLOT MEAN--###############
zlims=c(0,10)
pdf('d:/dropbox/working/diazotrophs/plots/nifh_maps.pdf',height=6,width=9)
par(mfrow=c(2,2),mar=c(2,2,2,5))
image(x=lon,y=lat,t(log10(rmean+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(tmean+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Tirchodesmium')
image(x=lon,y=lat,t(log10(uamean))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(ubmean+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-B')
dev.off()

# par(mfrow=c(2,2),mar=c(2,2,2,5))
# image(x=lon,y=lat,t(log10(rsd+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=c(0,9)); 
	# box()
	# map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	# image.plot(log10(rsd+1),legend.only=TRUE,col=viridis(20))
	# mtext('Richelia')
# image(x=lon,y=lat,t(log10(tsd+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=c(0,9)); 
	# box()
	# map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	# image.plot(log10(tsd+1),legend.only=TRUE,col=viridis(20))
	# mtext('Tirchodesmium')
# image(x=lon,y=lat,t(log10(uasd))[,90:1],xlab='',ylab='',col=viridis(20),zlim=c(0,9)); 
	# box()
	# map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	# image.plot(log10(uasd+1),legend.only=TRUE,col=viridis(20))
	# mtext('UCYN-A')
# image(x=lon,y=lat,t(log10(ubsd+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=c(0,9)); 
	# box()
	# map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	# image.plot(log10(ubmean),legend.only=TRUE,col=viridis(20))
	# mtext('UCYN-B')

##--PLOT nifH RANGE--######################
zlims <- c(0,10)
par(mfrow=c(2,2),mar=c(2,2,2,5))
image(x=lon,y=lat,t(log10(rrng+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(trng+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Tirchodesmium')
image(x=lon,y=lat,t(log10(uarng))[,90:1],xlab='',ylab='',col=viridis(20),zlim=c(0,10)); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(ubrng+1))[,90:1],xlab='',ylab='',col=viridis(20),zlim=c(0,10)); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-B')


#############################################################
## BIOMASS ANALYSIS #########################################
#############################################################
tlw  <- 1.22E-11; tup  <- 3.47E-8;  tcen <- (tlw+tup)/2
rlw  <- 5.75E-13; rup  <- 1.09E-11; rcen <- (rlw+rup)/2
ualw <- 6.59E-13; uaup <- 3.89E-9;  uacen <- (ualw+uaup)/2
ublw <- 1.39E-11; ubup <- 3.24E-10; ubcen <- (ublw+ubup)/2

Rcen <- rmean*rcen
Tcen <- tmean*tcen
UAcen <- uamean*uacen
UBcen <- ubmean*ubcen

Rrng  <- (rup*rmean - rlw*rmean)
Trng  <- (tup*tmean - tlw*tmean)
UArng <- (uaup*uamean - ualw*uamean)
UBrng <- (ubup*ubmean - ublw*ubmean)

##--PLOT MEAN--######################
pdf('d:/dropbox/working/diazotrophs/plots/biomass_maps.pdf',height=6,width=9)
zlims <- c(-10,2)
par(mfrow=c(2,2),mar=c(2,2,2,5))
image(x=lon,y=lat,t(log10(Rcen))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(Tcen))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Trichodesmium')
image(x=lon,y=lat,t(log10(uamean*uacen))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(ubmean*ubcen))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Trichodesmium')
dev.off()


##--PLOT RANGE--######################
zlims <- c(-10,2)
pdf('d:/dropbox/working/diazotrophs/plots/biomass_range_maps.pdf',height=6,width=9)
par(mfrow=c(2,2),mar=c(2,2,2,5))
image(x=lon,y=lat,t(log10(Rrng))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(Trng))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Trichodesmium')
image(x=lon,y=lat,t(log10(UArng))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(UBrng))[,90:1],xlab='',ylab='',col=viridis(20), zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Trichodesmium')
dev.off()

####################################################
## TWO COLUMN ######################################
####################################################


#####################################################################
## PACIFIC CENTERED PLOTS ###########################################
#####################################################################
loni <- c(70:180,1:69)

mapp <- function(lon,lat,map){
	#image(x=lon[1:90]+360,y=lat,t(map)[1:90,90:1],xlab='',ylab='',add=TRUE)
	#image(x=lon[101:180],y=lat,t(map)[101:180,90:1],xlab='',ylab='',add=TRUE)
	#image(x=lon[91:100]+360,y=lat,t(map)[91:100,90:1],add=TRUE)
	#image(x=lon+180,y=lat,t(map)[,90:1],xlab='',ylab='',add=TRUE)
	image(x=lon[1:90]+360,y=lat,t(map)[1:90,90:1],xlab='',ylab='',add=TRUE)
	#image(x=lon[91:180],y=lat,t(map)[91:180,90:1],xlab='',ylab='',add=TRUE)
	box(); axis(1);axis(2)
}

par(mfrow=c(2,2),mar=c(2,2,2,4),oma=c(2,2,2,8))
map('world2Hires', fill=TRUE, border=NA,col='grey',wrap=c(0,360),resolution=100)
	mapp(lon,lat,rmean); 
map('world2Hires', fill=TRUE, border=NA,col='grey',wrap=c(0,360)+20,resolution=100); 
	mapp(lon,lat,tmean)
map('world2Hires', fill=TRUE, border=NA,col='grey',wrap=c(0,360)+20,resolution=100); box(); axis(1);axis(2)
	mapp(lon,lat,uamean)
map('world2Hires', fill=TRUE, border=NA,col='grey',wrap=c(0,360)+20,resolution=100); box(); axis(1);axis(2)
	mapp(lon,lat,ubmean)
	image.plot(tmean,legend.only=TRUE,legend.args=list(line=2))


