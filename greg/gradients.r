
#diaz <- read.csv('~/dropbox/data/diazotrophs/maredat_diazotroph.csv',sep=',')
diaz <- read.csv('~/dropbox/working/DIAZOTROPHS/data/correction_integrated_04_23_2021.csv',sep=',')

refs <- unique(diaz$ref)

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




i <- 8
dd <- diaz[diaz$SOURCE..Data==refs[i],]
dd <- dd[order(dd$LATITUDE),]
	plot(dd$LATITUDE, dd$Richelia.nifH.Gene..copies.L.1,type='l',log='y')
	lines(dd$LATITUDE, dd$UCYN.B.nifH.Gene..copies.L.1.,col='red')
	lines(dd$LATITUDE, as.numeric(dd$UCYN.A1.nifH.Gene..copies.L.1.),col='blue')
	lines(dd$LATITUDE,dd$Trichodesmium.nifH.Gene..copies.L.1.,col='dark green')

	# for(i in 1:length(refs)){
# 	dd <- diaz[diaz$SOURCE..Data==refs[i],]
# 	dd <- dd[order(dd$LATITUDE),]
# 	plot(dd$LATITUDE,dd$Richelia.nifH.Gene..x106.copies.m.2.,type='l',log='y')
# 	lines(dd$LATITUDE,dd$UCYN.B.nifH.Gene..x106.copies.m.2.,col='red')
# 	lines(dd$LATITUDE,dd$UCYN.A.nifH.Gene..x106.copies.m.2.,col='blue')
# 	lines(dd$LATITUDE,dd$Trichodesmium.nifH.Gene..x106.copies.m.2.,col='dark green')
# }
# 

# for(i in 1:length(refs)){
# 	dd <- diaz[diaz$SOURCE..Data==refs[i],]
# 	dd <- dd[order(dd$LATITUDE),]
# 	plot(-999,ylim=c(1E-1,1E9),log='y',xlim=c(0,20))
# 	lines(dd$LATITUDE,dd$Richelia.nifH.Gene..x106.copies.m.2.+1)
# 	lines(dd$LATITUDE,dd$UCYN.B.nifH.Gene..x106.copies.m.2.+1,col='red')
# 	lines(dd$LATITUDE,dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1,col='blue')
# 	lines(dd$LATITUDE,dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1,col='dark green')
# }

# 
# for(i in 1:length(refs)){
# 	dd <- diaz[diaz$SOURCE..Data==refs[i],]
# 	dd <- dd[order(dd$LATITUDE),]
# 	plot(-999,ylim=c(1E-3,1E9),log='y',xlim=c(0,20))
# 	lines(dd$LATITUDE,(1/r_nifcell)*(dd$Richelia.nifH.Gene..x106.copies.m.2.+1))
# 	lines(dd$LATITUDE,(1/ucynb_nifcell)*(dd$UCYN.B.nifH.Gene..x106.copies.m.2.+1),col='red')
# 	lines(dd$LATITUDE,(1/ucyna_nifcell)*(dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1),col='blue')
# 	lines(dd$LATITUDE,(1/t_nifcell)*(dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1),col='dark green')
# }
# 


# dd <- diaz[diaz$SOURCE..Data==refs[31],]
# dd <- dd[order(dd$LATITUDE),]
dd <- diaz[diaz$ref==refs[31],]
dd <- dd[order(dd$lat),]

# x1 <- dd$Richelia.nifH.Gene..x106.copies.m.2.+1
# x2 <- dd$UCYN.B.nifH.Gene..x106.copies.m.2.+1
# x3 <- dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1
# x4 <- dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1

# x1 <- dd$Richelia.nifH.Gene..copies.L.1.+1
# x2 <- dd$UCYN.B.nifH.Gene..copies.L.1.+1
# x3 <- as.numeric(dd$UCYN.A1.nifH.Gene..copies.L.1.) +1
# x4 <- dd$Trichodesmium.nifH.Gene..copies.L.1.+1

x1 <- dd$richelia+1
x2 <- dd$ucyn_b+1
x3 <- dd$ucyn_a+1
x4 <- dd$trichodesmium+1



xx1 <- (1/r_nifcell)*x1
xx2 <- (1/ucynb_nifcell)*x2
xx3 <- (1/ucyna_nifcell)*x3
xx4 <- (1/t_nifcell)*x4
xx1_l <- (1/r_nifcell_l)*x1
xx2_l <- (1/ucynb_nifcell_l)*x2
xx3_l <- (1/ucyna_nifcell_l)*x3
xx4_l <- (1/t_nifcell_l)*x4
xx1_h <- (1/r_nifcell_h)*x1
xx2_h <- (1/ucynb_nifcell_h)*x2
xx3_h <- (1/ucyna_nifcell_h)*x3
xx4_h <- (1/t_nifcell_h)*x4

xxx1 <- r_Ccell*xx1
xxx2 <- ucynb_Ccell*xx2
xxx3 <- ucyna_Ccell*xx3
xxx4 <- t_Ccell*xx4
xxx1_l <- r_Ccell_l*xx1_h
xxx2_l <- ucynb_Ccell_l*xx2_h
xxx3_l <- ucyna_Ccell_l*xx3_h
xxx4_l <- t_Ccell_l*xx4_h
xxx1_h <- r_Ccell_h*xx1_l
xxx2_h <- ucynb_Ccell_h*xx2_l
xxx3_h <- ucyna_Ccell_h*xx3_l
xxx4_h <- t_Ccell_h*xx4_l

par(mfrow=c(1,1))
	plot(-999,ylim=c(1E-3,1E9),log='y',xlim=c(0,20))
	lines(dd$LATITUDE,x1+x2+x3+x4)
	lines(dd$LATITUDE,xx1+xx2+xx3+xx4)
	lines(dd$LATITUDE,xxx1+xxx2+xxx3+xxx4)


sumnif <- x1+x2+x3+x4
sumcell <- xx1+xx2+xx3+xx4
sumbio <- xxx1+xxx2+xxx3+xxx4

pdf('~/dropbox/working/diazotrophs/plots/diazotroph_gradients_04_16_2021.pdf',height=6,width=10)
#par(mfrow=c(2,3),mar=c(2,2,2,2),oma=c(2,2,2,5),xpd=TRUE)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(3,3,3,5),xpd=TRUE)
#   plot(-999,ylim=c(1E-2,1E9),log='y',xlim=c(0,20),las=1,ylab='')
# 	lines(dd$LATITUDE,x1+x2+x3+x4,lwd=2)
# 	lines(dd$LATITUDE,dd$Richelia.nifH.Gene..x106.copies.m.2.+1,col='orange',lty=2)
# 	lines(dd$LATITUDE,dd$UCYN.B.nifH.Gene..x106.copies.m.2.,col='red',lty=2)
# 	lines(dd$LATITUDE,dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1,col='blue',lty=2)
# 	lines(dd$LATITUDE,dd$Trichodesmium.nifH.Gene..x106.copies.m.2.,col='dark green',lty=2)
# 		mtext(side=2,expression('log(nifH [copies/m'^2*'])'),line=3)

	plot(-999,ylim=c(1E-3,1E9),log='y',xlim=c(0,20),las=1,ylab='')
	lines(dd$lat,xx1_h+xx2_h+xx3_h+xx4_h,lwd=2)
	lines(dd$lat,(1/r_nifcell_h)*(dd$richelia+1),col='orange',lty=2)
	lines(dd$lat,(1/ucynb_nifcell_h)*(dd$ucyn_b+1),col='red',lty=2)
	lines(dd$lat,(1/ucyna_nifcell_h)*(dd$ucyn_a+1),col='blue',lty=2)
	lines(dd$lat,(1/t_nifcell_h)*(dd$trichodesmium+1),col='dark green',lty=2)
		mtext(side=2,expression('log10(cells [cells/m'^2*'])'),line=3.5)
		legend('topleft',legend=c('Total','Richelia','UCYN A','UCYN B','Trichodesmium'),lty=c(1,2,2,2,2),lwd=c(2,1,1,1,1),
		       col=c('black','orange','blue','red','dark green'),bty='n',cex=0.7)
		
		plot(-999,ylim=c(1E-3,1E9),log='y',xlim=c(0,20),las=1,ylab='')
		lines(dd$lat,xx1_l+xx2_l+xx3_l+xx4_l,lwd=2)
		lines(dd$lat,(1/r_nifcell_l)*(dd$richelia+1),col='orange',lty=2)
		lines(dd$lat,(1/ucynb_nifcell_l)*(dd$ucyn_b+1),col='red',lty=2)
		lines(dd$lat,(1/ucyna_nifcell_l)*(dd$ucyn_a+1),col='blue',lty=2)
		lines(dd$lat,(1/t_nifcell_l)*(dd$trichodesmium+1),col='dark green',lty=2)
		
		
	plot(-999,ylim=c(1E-13,5E0),log='y',xlim=c(0,20),las=1,ylab='')
	lines(dd$lat,xxx1_l+xxx2_l+xxx3_l+xxx4_l,lwd=2)
	lines(dd$lat,(1/r_nifcell_h)*(r_Ccell_l)*(dd$richelia+1),col='orange',lty=2)
	lines(dd$lat,(1/ucynb_nifcell_h)*(ucynb_Ccell_l)*(dd$ucyn_b+1),col='red',lty=2)
	lines(dd$lat,(1/ucyna_nifcell_h)*(ucyna_Ccell_l)*(dd$ucyn_a+1),col='blue',lty=2)
	lines(dd$lat,(1/t_nifcell_h)*(t_Ccell_l)*(dd$trichodesmium+1),col='dark green',lty=2)
		mtext(side=2,expression('log10(carbon [mmol/m'^2*'])'),line=3.5)


# plot(-999,ylim=c(1E-2,1E9),log='y',xlim=c(0,20),las=1,ylab='')
# lines(dd$LATITUDE,x1+x2+x3+x4,lwd=2)
# lines(dd$LATITUDE,dd$Richelia.nifH.Gene..x106.copies.m.2.+1,col='orange',lty=2)
# lines(dd$LATITUDE,dd$UCYN.B.nifH.Gene..x106.copies.m.2.,col='red',lty=2)
# lines(dd$LATITUDE,dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1,col='blue',lty=2)
# lines(dd$LATITUDE,dd$Trichodesmium.nifH.Gene..x106.copies.m.2.,col='dark green',lty=2)


plot(-999,ylim=c(1E-13,5E0),log='y',xlim=c(0,20),las=1,ylab='')
lines(dd$lat,xxx1_h+xxx2_h+xxx3_h+xxx4_h,lwd=2)
lines(dd$lat,(1/r_nifcell_l)*(r_Ccell_h)*(dd$richelia+1),col='orange',lty=2)
lines(dd$lat,(1/ucynb_nifcell_l)*(ucynb_Ccell_h)*(dd$ucyn_b+1),col='red',lty=2)
lines(dd$lat,(1/ucyna_nifcell_l)*(ucyna_Ccell_h)*(dd$ucyn_a+1),col='blue',lty=2)
lines(dd$lat,(1/t_nifcell_l)*(t_Ccell_h)*(dd$trichodesmium+1),col='dark green',lty=2)
mtext(side=1,outer=TRUE,expression('Latitude ('*degree*'N)'),line=0.5)

dev.off()

############################################################
A <- data.frame(x=rnorm(100, 20, 2), y=rnorm(100, 20, 2))
B <- data.frame(x=rnorm(100, 21, 1), y=rnorm(100, 21, 1))

# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

# Plot both groups
plot(y ~ x, A, ylim=range(c(A$y, B$y)), xlim=range(c(A$x, B$x)), pch=1,
               main="Scatter plot of two groups")
points(y ~ x, B, pch=3)

# Add legend to top right, outside plot region
legend("topright", inset=c(-0.2,0), legend=c("A","B"), pch=c(1,3), title="Group")



ddd <- dd[order(dd$LONGITUDE),]

x1 <- ddd$Richelia.nifH.Gene..x106.copies.m.2.+1
x2 <- ddd$UCYN.B.nifH.Gene..x106.copies.m.2.+1
x3 <- ddd$UCYN.A.nifH.Gene..x106.copies.m.2.+1
x4 <- ddd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1

xx1 <- (1/r_nifcell)*x1
xx2 <- (1/ucynb_nifcell)*x2
xx3 <- (1/ucyna_nifcell)*x3
xx4 <- (1/t_nifcell)*x4

xxx1 <- r_Ccell*xx1
xxx2 <- ucynb_Ccell*xx2
xxx3 <- ucyna_Ccell*xx3
xxx4 <- t_Ccell*xx4

par(mfrow=c(1,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
	plot(-999,ylim=c(1E-2,1E9),log='y',xlim=c(-70,-20),las=1,ylab='')
	lines(ddd$LONGITUDE,x1+x2+x3+x4,lwd=2)
	lines(ddd$LONGITUDE,ddd$Richelia.nifH.Gene..x106.copies.m.2.+1,col='orange',lty=2)
	lines(ddd$LONGITUDE,ddd$UCYN.B.nifH.Gene..x106.copies.m.2.,col='red',lty=2)
	lines(ddd$LONGITUDE,ddd$UCYN.A.nifH.Gene..x106.copies.m.2.+1,col='blue',lty=2)
	lines(ddd$LONGITUDE,ddd$Trichodesmium.nifH.Gene..x106.copies.m.2.,col='dark green',lty=2)
		mtext(adj=0.1,expression('log10(nifH [copies/m'^2*'])'))

	plot(-999,ylim=c(1E-3,1E8),log='y',xlim=c(-70,-20),las=1)
	lines(ddd$LONGITUDE,xx1+xx2+xx3+xx4,lwd=2)
	lines(ddd$LONGITUDE,(1/r_nifcell)*(ddd$Richelia.nifH.Gene..x106.copies.m.2.+1),col='orange',lty=2)
	lines(ddd$LONGITUDE,(1/ucynb_nifcell)*(ddd$UCYN.B.nifH.Gene..x106.copies.m.2.+1),col='red',lty=2)
	lines(ddd$LONGITUDE,(1/ucyna_nifcell)*(ddd$UCYN.A.nifH.Gene..x106.copies.m.2.+1),col='blue',lty=2)
	lines(dd$LATITUDE,(1/t_nifcell)*(dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1),col='dark green',lty=2)
		mtext(adj=0.1,expression('log10(cells [cells/m'^2*'])'))

	plot(-999,ylim=c(1E-13,1E9),log='y',xlim=c(-70,-20),las=1)
	lines(ddd$LONGITUDE,xxx1+xxx2+xxx3+xxx4,lwd=2)
	lines(dd$LATITUDE,(1/r_nifcell)*(r_Ccell)*(dd$Richelia.nifH.Gene..x106.copies.m.2.+1),col='orange',lty=2)
	lines(dd$LATITUDE,(1/ucynb_nifcell)*(ucynb_Ccell)*(dd$UCYN.B.nifH.Gene..x106.copies.m.2.+1),col='red',lty=2)
	lines(dd$LATITUDE,(1/ucyna_nifcell)*(ucyna_Ccell)*(dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1),col='blue',lty=2)
	lines(dd$LATITUDE,(1/t_nifcell)*(t_Ccell)*(dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1),col='dark green',lty=2)
		mtext(adj=0.1,expression('log10(carbon [mmol/m'^2*'])'))
legend('topright',legend=c('Total','Richelia','UCYN A','UCYN B','Trichodesmium'),lty=c(1,2,2,2,2),lwd=c(2,1,1,1,1),
	col=c('black','orange','blue','red','dark green'),bty='n')
mtext('Latitude',outer=TRUE,side=1,line=0.5)












par(mfrow=c(1,3))
plot(-1e14,ylim=c(1E-3,1E8),xlim=c(0,20),las=1,ylab='')
	lines(dd$LATITUDE,dd$Richelia.nifH.Gene..x106.copies.m.2.)
	lines(dd$LATITUDE,dd$UCYN.B.nifH.Gene..x106.copies.m.2.,col='red')
	lines(dd$LATITUDE,dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1,col='blue')
	lines(dd$LATITUDE,dd$Trichodesmium.nifH.Gene..x106.copies.m.2.,col='dark green')

plot(-999,ylim=c(1E-3,1E8),xlim=c(0,20))
	lines(dd$LATITUDE,(1/r_nifcell)*(dd$Richelia.nifH.Gene..x106.copies.m.2.+1))
	lines(dd$LATITUDE,(1/ucynb_nifcell)*(dd$UCYN.B.nifH.Gene..x106.copies.m.2.+1),col='red')
	lines(dd$LATITUDE,(1/ucyna_nifcell)*(dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1),col='blue')
	lines(dd$LATITUDE,(1/t_nifcell)*(dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1),col='dark green')

	plot(-999,ylim=c(1E-13,1E9),log='y',xlim=c(0,20))
	lines(dd$LATITUDE,(1/r_nifcell)*(r_Ccell)*(dd$Richelia.nifH.Gene..x106.copies.m.2.+1))
	lines(dd$LATITUDE,(1/ucynb_nifcell)*(ucynb_Ccell)*(dd$UCYN.B.nifH.Gene..x106.copies.m.2.+1),col='red')
	lines(dd$LATITUDE,(1/ucyna_nifcell)*(ucyna_Ccell)*(dd$UCYN.A.nifH.Gene..x106.copies.m.2.+1),col='blue')
	lines(dd$LATITUDE,(1/t_nifcell)*(t_Ccell)*(dd$Trichodesmium.nifH.Gene..x106.copies.m.2.+1),col='dark green')



