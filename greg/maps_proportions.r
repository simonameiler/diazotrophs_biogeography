library(R.matlab)
library(fields)
library(maps)
library(mapdata)
library(viridis)
data(world2HiresMapEnv)
data(world2LoresMapEnv)

######################################################
## READ IN DATA ######################################
######################################################
setwd('~/dropbox/working/diazotrophs/tang-cassar-2019/')

D <- readMat('Diazotrophs_RF_global_monthly.mat')  #read the full dataset

r  <- D$Richelia.RF.global.monthly  #extract individual groups
t  <- D$Tricho.RF.global.monthly
ua <- D$UCYNA.RF.global.monthly
ub <- D$UCYNB.RF.global.monthly
tot <- r + t + ua + ub

r[r=='1'] <- 0    #set 1s to true zeros
t[t=='1'] <- 0
ua[ua=='1'] <- 0
ub[ub=='1'] <- 0

lon <- seq(-180,180,length.out=180)
lat <- seq(-90,90,length.out=90)

rmean  <- apply(r,c(1,2),mean)
tmean  <- apply(t,c(1,2),mean)
uamean <- apply(ua,c(1,2),mean)
ubmean <- apply(ub,c(1,2),mean)

rrng   <- apply(r,c(1,2),function(x) diff(range(x)))
trng   <- apply(t,c(1,2),function(x) diff(range(x)))
uarng  <- apply(ua,c(1,2),function(x) diff(range(x)))
ubrng  <- apply(ub,c(1,2),function(x) diff(range(x)))
totrng <- apply(tot,c(1,2),function(x) diff(range(x)))

#########################################################################
## CONVERSIONS ##########################################################
#########################################################################
t_nifcell_h <- 500
t_nifcell_l <- 100
t_nifcell   <- mean(c(t_nifcell_h,t_nifcell_l))

ucyna_nifcell_h <- 27.91
ucyna_nifcell_l <- 0.49
ucyna_nifcell   <- mean(c(ucyna_nifcell_h,ucyna_nifcell_l))

ucynb_nifcell_h <- 3.6
ucynb_nifcell_l <- 3.6
ucynb_nifcell   <- mean(c(ucynb_nifcell_h,ucynb_nifcell_l))

r_nifcell_h <- 1000
r_nifcell_l <- 100
r_nifcell   <- mean(c(r_nifcell_h,r_nifcell_l))


t_Ccell_h <- 3.47E-6
t_Ccell_l <- 6.08E-9
t_Ccell   <- mean(t_Ccell_h,t_Ccell_l) 

ucyna_Ccell_h <- 1.91E-9
ucyna_Ccell_l <- 1.83E-11
ucyna_Ccell   <- mean(ucyna_Ccell_h,ucyna_Ccell_l)

ucynb_Ccell_h <- 1.16E-9
ucynb_Ccell_l <- 5E-11
ucynb_Ccell   <- mean(ucynb_Ccell_h,ucynb_Ccell_l)

r_Ccell_h <- 1.09E-9
r_Ccell_l <- 5.75E-10
r_Ccell   <- mean(r_Ccell_h,r_Ccell_l)

###########################################################################
## MEAN CELLULAR ABUNDANCE BY SPECIES #####################################
###########################################################################
rcell <- t(rmean*(1/r_nifcell))[,90:1]
tcell <- t(tmean*(1/t_nifcell))[,90:1]
uacell <- t(uamean*(1/ucyna_nifcell))[,90:1]
ubcell <- t(ubmean*(1/ucynb_nifcell))[,90:1]

rcell_l <- t(rmean*(1/r_nifcell_l))[,90:1]
tcell_l <- t(tmean*(1/t_nifcell_l))[,90:1]
uacell_l <- t(uamean*(1/ucyna_nifcell_l))[,90:1]
ubcell_l <- t(ubmean*(1/ucynb_nifcell_l))[,90:1]

rcell_h <- t(rmean*(1/r_nifcell_h))[,90:1]
tcell_h <- t(tmean*(1/t_nifcell_h))[,90:1]
uacell_h <- t(uamean*(1/ucyna_nifcell_h))[,90:1]
ubcell_h <- t(ubmean*(1/ucynb_nifcell_h))[,90:1]

rcelltot <- r*(1/r_nifcell)
tcelltot <- t*(1/t_nifcell)
uacelltot <- ua*(1/ucyna_nifcell)
ubcelltot <- ub*(1/ucynb_nifcell)

############################################################
## PLOT MEAN CELL ABUNDANCE BY SPECIES #####################
############################################################
cellr  <- t(rmean*(1/r_nifcell))[,90:1]
cellt  <- t(tmean*(1/t_nifcell))[,90:1]
cellua <- t(uamean*(1/ucyna_nifcell))[,90:1]
cellub <- t(ubmean*(1/ucynb_nifcell))[,90:1]

celltot <- cellr + cellt + cellua + cellub

rC  <- rcell*r_Ccell
tC  <- tcell*t_Ccell
uaC <- uacell*ucyna_Ccell
ubC <- ubcell*ucynb_Ccell

Ctot <- rC + tC + uaC + ubC


zlims=c(0,1) 
#pdf('~/dropbox/working/diazotrophs/plots/mean_cell_04_01_2020.pdf',height=9,width=8.5)
CairoPDF('~/dropbox/working/diazotrophs/plots/mean_cell_04_01_2020.pdf',height=9,width=8.5)
par(mfrow=c(4,2),mar=c(2,2,2,5),cex.axis=0.8,oma=c(2,2,3,2),xpd=FALSE)
image(x=lon,y=lat,cellr/celltot ,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
	mtext('Proportion of Cellular Concentration',line=2)
image(x=lon,y=lat,rC/Ctot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
	mtext('Proportion of Carbon Concentration',line=2)
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))

image(x=lon,y=lat,cellt/celltot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('Trichodesmium')
image(x=lon,y=lat,tC/Ctot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('Trichodesmium')
image(x=lon,y=lat,cellua/celltot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('UCYN-A')
image(x=lon,y=lat,uaC/Ctot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('UCYN-A')

image(x=lon,y=lat,cellub/celltot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('UCYN-B')
image(x=lon,y=lat,ubC/Ctot,xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('UCYN-B')

dev.off()
