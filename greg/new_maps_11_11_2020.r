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
setwd('d:/dropbox/working/diazotrophs/tang-cassar-2019/')

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
########################################################################
## PLOT MEAN NIFH ######################################################
########################################################################
zlims=c(0,10)
par(mfrow=c(2,2),mar=c(2,2,2,5))
image(x=lon,y=lat,t(log10(rmean))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(tmean))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Tirchodesmium')
image(x=lon,y=lat,t(log10(uamean))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(ubmean))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-B')

zlims=c(0,10)
par(mfrow=c(2,2),mar=c(2,2,2,5))
image(x=lon,y=lat,t(log10(rrng))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(trng))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Tirchodesmium')
image(x=lon,y=lat,t(log10(uarng))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(ubrng))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-B')

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

r_Ccell_h <- 8.16E-7
r_Ccell_l <- 1.52E-9
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
zlims=c(-3,9) 
pdf('d:/dropbox/working/diazotrophs/plots/mean_cell_12_01_2020.pdf',height=5,width=8.5)
par(mfrow=c(2,2),mar=c(2,2,2,5),cex.axis=0.8,oma=c(2,2,3,2))
image(x=lon,y=lat,t(log10(rmean*(1/r_nifcell)))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,t(log10(tmean*(1/t_nifcell)))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Tirchodesmium')
image(x=lon,y=lat,t(log10(uamean*(1/ucyna_nifcell)))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('UCYN-A')
image(x=lon,y=lat,t(log10(ubmean*(1/ucynb_nifcell)))[,90:1],xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-B')
	mtext(outer=TRUE,'Cell Abundance')
dev.off()

############################################################
## PLOT MEAN CARBON BY SPECIES #####################
############################################################
rC  <- rcell*r_Ccell
tC  <- tcell*t_Ccell
uaC <- uacell*ucyna_Ccell
ubC <- ubcell*ucynb_Ccell

zlims=c(-10,2) 
pdf('d:/dropbox/working/diazotrophs/plots/mean_carbon_12_01_2020.pdf',height=5,width=8.5)
par(mfrow=c(2,2),mar=c(2,2,2,5),cex.axis=0.8,oma=c(2,2,3,2))
image(x=lon,y=lat,log10(rC),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Richelia')
image(x=lon,y=lat,log10(tC),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('Trichodesmium')
image(x=lon,y=lat,log10(uaC),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	mtext('UCYN-A')
image(x=lon,y=lat,log10(ubC),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('UCYN-B')
	mtext(outer=TRUE,'Carbon')
dev.off()







###############################################################
## TOTAL CELL + CARBON ABUNDANCE ##############################
###############################################################
cell <- rcell + tcell + uacell + ubcell
 
cell_l <- t(rmean*(1/r_nifcell_l))[,90:1] + t(tmean*(1/t_nifcell_l))[,90:1] + 
	t(uamean*(1/ucyna_nifcell_l))[,90:1] + t(ubmean*(1/ucynb_nifcell_l))[,90:1]
cell_h <- t(rmean*(1/r_nifcell_h))[,90:1] + t(tmean*(1/t_nifcell_h))[,90:1] + 
	t(uamean*(1/ucyna_nifcell_h))[,90:1] + t(ubmean*(1/ucynb_nifcell_h))[,90:1]


rCtot  <- rcelltot*r_Ccell
tCtot  <- tcelltot*t_Ccell
uaCtot <- uacelltot*ucyna_Ccell
ubCtot <- ubcelltot*ucynb_Ccell


C   <- rC + tC + uaC + ubC 
C_l <- rcell*r_Ccell_l + tcell*t_Ccell_l + uacell*ucyna_Ccell_l + ubcell*ucynb_Ccell_l
C_h <- rcell*r_Ccell_h + tcell*t_Ccell_h + uacell*ucyna_Ccell_h + ubcell*ucynb_Ccell_h

C_l_l <- rcell_h*r_Ccell_l + tcell_h*t_Ccell_l + uacell_h*ucyna_Ccell_l + ubcell_h*ucynb_Ccell_l
C_h_h <- rcell_l*r_Ccell_h + tcell_l*t_Ccell_h + uacell_l*ucyna_Ccell_h + ubcell_l*ucynb_Ccell_h


#######################################################################
## LATITUDINAL PLOTS ##################################################
#######################################################################

pdf('d:/dropbox/working/diazotrophs/plots/latitudinal.pdf',height=4,width=7)
lats <- seq(-90,90,length.out=90)
ylims=c(-10,0)
par(mfrow=c(2,2),mar=c(1,1,2,2),oma=c(3,3,2,2))
plot(lats,colMeans(log10(rC),na.rm=TRUE),type='l',ylim=ylims)
	lines(lats,colMeans(log10(rcell_h*r_Ccell_l),na.rm=TRUE),lty=2)
	lines(lats,colMeans(log10(rcell_l*r_Ccell_h),na.rm=TRUE),lty=2)
	mtext('Richelia')
	
plot(lats,colMeans(log10(tC),na.rm=TRUE),type='l',ylim=ylims)
	lines(lats,colMeans(log10(tcell_h*t_Ccell_l),na.rm=TRUE),type='l',lty=2)
	lines(lats,colMeans(log10(tcell_l*t_Ccell_h)+0.1,na.rm=TRUE),type='l',lty=2)
	mtext('Trichodesmium')

plot(lats,colMeans(log10(uaC),na.rm=TRUE),type='l',ylim=ylims)
	lines(lats,colMeans(log10(uacell_h*ucyna_Ccell_l),na.rm=TRUE),type='l',lty=2)
	lines(lats,colMeans(log10(uacell_l*ucyna_Ccell_h)-0.2,na.rm=TRUE),type='l',lty=2)
	mtext('UCYN-A')

plot(colMeans(log10(ubC),na.rm=TRUE),type='l',ylim=ylims)
	lines(colMeans(log10(ubcell_h*ucynb_Ccell_l),na.rm=TRUE),type='l',lty=2)
	lines(colMeans(log10(ubcell_l*ucynb_Ccell_h)+0.8,na.rm=TRUE),type='l',lty=2)
	mtext('UCYN-B')
mtext(side=1,outer=TRUE,'Latitude',line=1)
mtext(side=2,outer=TRUE,expression('log10(Carbon Concentration [mgC/m'^3*']'),line=1)
mtext(outer=TRUE,'Carbon')
dev.off()


##################################################################
##--TOTAL CELL AND TOTAL CARBON 2 x 3 PANEL ######################
##################################################################
zlims=c(-0,10)
pdf('d:/dropbox/working/diazotrophs/plots/total_cell_total_carbon_2x3_12_01_2020.pdf',height=5,width=12)
par(mfrow=c(2,3),mar=c(2,2,2,5))
image(x=lon,y=lat,log10(cell),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('log10(Cell Concentration)')
image(x=lon,y=lat,log10(cell_l-cell_h),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('log10(Cell Concentration Range)')
#image(x=lon,y=lat,(cell_l-cell_h)/cell,xlab='',ylab='',col=viridis(20)); 
#	box()
#	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
#	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
#	#mtext('UCYN-B')

plot(-999,type='n',xaxt='n',yaxt='n',bty='n')

zlims=c(-9,2)
image(x=lon,y=lat,log10(C),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('log10(Carbon Concentration)')
image(x=lon,y=lat,log10(C_h-C_l),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('log10(Carbon Concentration Range)')

image(x=lon,y=lat,log10(C_h_h-C_l_l),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	mtext('log10(Carbon Concentration Range)')
dev.off()

#############################################################
C   <- rC + tC + uaC + ubC 
C_l <- rcell*r_Ccell_l + tcell*t_Ccell_l + uacell*ucyna_Ccell_l + ubcell*ucynb_Ccell_l
C_h <- rcell*r_Ccell_h + tcell*t_Ccell_h + uacell*ucyna_Ccell_h + ubcell*ucynb_Ccell_h

C_l_l <- rcell_h*r_Ccell_l + tcell_h*t_Ccell_l + uacell_h*ucyna_Ccell_l + ubcell_h*ucynb_Ccell_l
C_h_h <- rcell_l*r_Ccell_h + tcell_l*t_Ccell_h + uacell_l*ucyna_Ccell_h + ubcell_l*ucynb_Ccell_h

#####################################################################
## MEAN CARBON SEASONAL RANGE #######################################
#####################################################################
rCrng <- t(rrng)[,90:1]*(1/r_nifcell_l)*r_Ccell
tCrng <- t(trng)[,90:1]*(1/t_nifcell_l)*t_Ccell
uaCrng <- t(uarng)[,90:1]*(1/ucyna_nifcell_l)*ucyna_Ccell
ubCrng <- t(ubrng)[,90:1]*(1/ucynb_nifcell_l)*ucynb_Ccell


image(x=lon,y=lat,log10(tCrng),xlab='',ylab='',col=viridis(20)); 
#image(x=lon,y=lat,(cell_l-cell_h)/celltotrng,xlab='',ylab='',col=viridis(20),zlim=c(0,1)); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	#mtext('UCYN-B')

C_l <- rcell*r_Ccell_l + tcell*t_Ccell_l + uacell*ucyna_Ccell_l + ubcell*ucynb_Ccell_l
C_h <- rcell*r_Ccell_h + tcell*t_Ccell_h + uacell*ucyna_Ccell_h + ubcell*ucynb_Ccell_h

C_l_l <- rcell_h*r_Ccell_l + tcell_h*t_Ccell_l + uacell_h*ucyna_Ccell_l + ubcell_h*ucynb_Ccell_l
C_h_h <- rcell_l*r_Ccell_h + tcell_l*t_Ccell_h + uacell_l*ucyna_Ccell_h + ubcell_l*ucynb_Ccell_h

###################################################################
## SPECIES CARBON UNCERTAINTY DIVIDED BY SEASONAL RANGE ###########
###################################################################
zlims <- c(0,1.2) 
pdf('d:/dropbox/working/diazotrophs/plots/meanspeciesCuncertainy_seasonal.pdf',height=5,width=8)
par(mfrow=c(2,2),mar=c(2,2,2,5))
rC_hl <- rcell_l*r_Ccell_h - rcell_h*r_Ccell_l 
	image(x=lon,y=lat,log10(rCrng/rC_hl),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
tC_hl <- tcell_l*t_Ccell_h - tcell_h*t_Ccell_l 
	image(x=lon,y=lat,log10(tCrng/tC_hl),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
uaC_hl <- uacell_l*ucyna_Ccell_h - uacell_h*ucyna_Ccell_l 
	image(x=lon,y=lat,log10(uaCrng/uaC_hl),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
ubC_hl <- ubcell_l*ucynb_Ccell_h - ubcell_h*ucynb_Ccell_l 
	image(x=lon,y=lat,log10(ubCrng/ubC_hl),xlab='',ylab='',col=viridis(20),zlim=zlims); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
dev.off()

celltot <- rcelltot + tcelltot + uacelltot + ubcelltot
celltotrng <- t(apply(celltot,c(1,2),function(x) diff(range(x))))[,90:1]

Ctot <- rCtot + tCtot + uaCtot + ubCtot
Ctotrng <- t(apply(Ctot,c(1,2),function(x) diff(range(x))))[,90:1]


####################################################################
## TOTAL CELL UNDERTAINTY DIVIDED BY SEASONAL RANGE ################
####################################################################
image(x=lon,y=lat,log10((cell_l-cell_h)/celltotrng),xlab='',ylab='',col=viridis(20),zlim=c(-3,2)); 
#image(x=lon,y=lat,(cell_l-cell_h)/celltotrng,xlab='',ylab='',col=viridis(20),zlim=c(0,1)); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	#mtext('UCYN-B')
image(x=lon,y=lat,log10((C_h-C_l)/Ctotrng),xlab='',ylab='',col=viridis(20),zlim=c(-1.5,0.5)); 
#image(x=lon,y=lat,(cell_l-cell_h)/celltotrng,xlab='',ylab='',col=viridis(20),zlim=c(0,1)); 
	box()
	map(add=TRUE,fill=TRUE,resolution=1000,col='grey')
	#image.plot(matrix(zlims),legend.only=TRUE,col=viridis(20))
	#mtext('UCYN-B')

	
################################################	
