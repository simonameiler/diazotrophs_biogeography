
#diaz <- read.csv('~/dropbox/data/diazotrophs/maredat_diazotroph.csv',sep=',')
#diaz <- read.csv('~/dropbox/working/DIAZOTROPHS/data/correction_integrated_04_23_2021.csv',sep=',')
diaz <- read.csv('~/dropbox/working/DIAZOTROPHS/data/correction_integrated_05_28_2021.csv',sep=',')
refs <- unique(diaz$ref)

############################################################
## CONVERSIONS #############################################
############################################################
t_nifcell_h <- 1012.45
t_nifcell_l <- 2.64
  t_nifcell   <- mean(c(t_nifcell_h,t_nifcell_l))
ucyna_nifcell_h <- 27.81
ucyna_nifcell_l <- 0.49
  ucyna_nifcell   <- mean(c(ucyna_nifcell_h,ucyna_nifcell_l))
ucynb_nifcell_h <- 3.6
ucynb_nifcell_l <- 3.6
  ucynb_nifcell   <- mean(c(ucynb_nifcell_h,ucynb_nifcell_l))
r_nifcell_h <- 401.56
r_nifcell_l <- 1.35
  r_nifcell   <- mean(c(r_nifcell_h,r_nifcell_l))
t_Ccell_h <- 4.75E-8
t_Ccell_l <- 1.48E-9
  t_Ccell   <- mean(t_Ccell_h,t_Ccell_l) 
ucyna_Ccell_h <- 1.91E-9
ucyna_Ccell_l <- 1.86E-11
  ucyna_Ccell   <- mean(ucyna_Ccell_h,ucyna_Ccell_l)
ucynb_Ccell_h <- 1.17E-9
ucynb_Ccell_l <- 5E-11
  ucynb_Ccell   <- mean(ucynb_Ccell_h,ucynb_Ccell_l)
r_Ccell_h <- 2.56E-9
r_Ccell_l <- 5.73E-11
  r_Ccell   <- mean(r_Ccell_h,r_Ccell_l)


dd <- diaz[diaz$ref==refs[31],]
dd <- dd[order(dd$lat),]

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

pdf('~/dropbox/working/diazotrophs/plots/diazotroph_gradients_05_28_2021.pdf',height=6,width=10)
#par(mfrow=c(2,3),mar=c(2,2,2,2),oma=c(2,2,2,5),xpd=TRUE)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(3,3,3,5),xpd=TRUE)

	plot(-999,ylim=c(1E-3,1E6),log='y',xlim=c(0,20),las=1,ylab='')
	lines(dd$lat,xx1_h+xx2_h+xx3_h+xx4_h,lwd=2)
	lines(dd$lat,(1/r_nifcell_h)*(dd$richelia+1),col='orange',lty=2)
	lines(dd$lat,(1/ucynb_nifcell_h)*(dd$ucyn_b+1),col='red',lty=2)
	lines(dd$lat,(1/ucyna_nifcell_h)*(dd$ucyn_a+1),col='blue',lty=2)
	lines(dd$lat,(1/t_nifcell_h)*(dd$trichodesmium+1),col='dark green',lty=2)
		mtext(side=2,expression('Cells [cells/m'^2*']'),line=3.5)
		legend('topleft',legend=c('Total','Richelia','UCYN A','UCYN B','Trichodesmium'),lty=c(1,2,2,2,2),lwd=c(2,1,1,1,1),
		       col=c('black','orange','blue','red','dark green'),bty='n',cex=0.7)
		mtext(adj=0,'a)')

		plot(-999,ylim=c(1E-3,1E6),log='y',xlim=c(0,20),las=1,ylab='')
		lines(dd$lat,xx1_l+xx2_l+xx3_l+xx4_l,lwd=2)
		lines(dd$lat,(1/r_nifcell_l)*(dd$richelia+1),col='orange',lty=2)
		lines(dd$lat,(1/ucynb_nifcell_l)*(dd$ucyn_b+1),col='red',lty=2)
		lines(dd$lat,(1/ucyna_nifcell_l)*(dd$ucyn_a+1),col='blue',lty=2)
		lines(dd$lat,(1/t_nifcell_l)*(dd$trichodesmium+1),col='dark green',lty=2)
		mtext(adj=0,'b)')
		
		
	plot(-999,ylim=c(1E-13,1E-3),log='y',xlim=c(0,20),las=1,ylab='')
	lines(dd$lat,xxx1_l+xxx2_l+xxx3_l+xxx4_l,lwd=2)
	lines(dd$lat,(1/r_nifcell_h)*(r_Ccell_l)*(dd$richelia+1),col='orange',lty=2)
	lines(dd$lat,(1/ucynb_nifcell_h)*(ucynb_Ccell_l)*(dd$ucyn_b+1),col='red',lty=2)
	lines(dd$lat,(1/ucyna_nifcell_h)*(ucyna_Ccell_l)*(dd$ucyn_a+1),col='blue',lty=2)
	lines(dd$lat,(1/t_nifcell_h)*(t_Ccell_l)*(dd$trichodesmium+1),col='dark green',lty=2)
		mtext(side=2,expression('Carbon [mmol/m'^2*']'),line=3.5)
		mtext(adj=0,'c)')

  plot(-999,ylim=c(1E-13,1E-3),log='y',xlim=c(0,20),las=1,ylab='')
  lines(dd$lat,xxx1_h+xxx2_h+xxx3_h+xxx4_h,lwd=2)
  lines(dd$lat,(1/r_nifcell_l)*(r_Ccell_h)*(dd$richelia+1),col='orange',lty=2)
  lines(dd$lat,(1/ucynb_nifcell_l)*(ucynb_Ccell_h)*(dd$ucyn_b+1),col='red',lty=2)
  lines(dd$lat,(1/ucyna_nifcell_l)*(ucyna_Ccell_h)*(dd$ucyn_a+1),col='blue',lty=2)
  lines(dd$lat,(1/t_nifcell_l)*(t_Ccell_h)*(dd$trichodesmium+1),col='dark green',lty=2)
  mtext(side=1,outer=TRUE,expression('Latitude ('*degree*'N)'),line=0.5)
  mtext(adj=0,'d)')

dev.off()

