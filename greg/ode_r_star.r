library(deSolve) #package for ODE solver
library(colorRamps)
library(RColorBrewer)
parameters <- list(N_0  =1.6,      #inflow N concentration
				P_0  =2*1.6/16, #inflow P concentration
				Fe_0 =1E-6,     #inflow Fe concentration
				f_atm=1E-1,     #atmospheric deposition flux
				kappa=0.1,      #mixing coefficient
				mu_p =2.5,      #max growth rate non-diazotrophs
				k_pN =5.6E-2,   #half-sat non-diaz wrt N	 
				k_pP =3.5E-2,   #half-sat non-diaz wrt P
				k_pFe=3.5E-4,   #
				r_pP =0.0625,
				r_pFe=6.25E-5,
				mu_d =1.25,
				k_dP =3.5E-2,
				k_dFe=1.1E-3,
				r_dP =0.025,
				r_dFe=7.5E-4,
				m    =0.05)       

state <- c(B_p = 1, #initial value non-diaz biomass
           B_d = 1, #initial value diaz biomass
           N   = 1, #initial value N
           P   = 1, #initial value P
		   Fe  = 1) #initial value Fe

dxdt <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
	gamma_p = min(N/(N+k_pN), P/(P+k_pP), Fe/(Fe+k_pFe))   #nutrient limitation for non-diaz
	gamma_d = min(            P/(P+k_dP), Fe/(Fe+k_dFe))   #nutrient limitation for diaz
	
	dB_p =  (mu_p*gamma_p - m - kappa)*B_p                 #incremental change in non-diaz
	dB_d =  (mu_d*gamma_d - m - kappa)*B_d                 #incremental change in diaz
	dN   = -mu_p*gamma_p*B_p                                + m*B_p       + m*B_d       + kappa*(N_0-N) #incremental change in N
	dP   = -mu_p*gamma_p*B_p*r_pP  - mu_d*gamma_d*B_d*r_dP  + m*B_p*r_pP  + m*B_d*r_dP  + kappa*(P_0-P)            #incremental change in P
	dFe  = -mu_p*gamma_p*B_p*r_pFe - mu_d*gamma_d*B_d*r_dFe + m*B_p*r_pFe + m*B_d*r_dFe + kappa*(Fe_0-Fe) + f_atm  #incremental change in Fe
    
    list(c(dB_p,dB_d,dN,dP,dFe))})}

dxdt_pref <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
	gamma_p = min(N/(N+k_pN), P/(P+k_pP), Fe/(Fe+k_pFe))   #nutrient limitation for non-diaz
	gamma_d = min(            P/(P+k_dP), Fe/(Fe+k_dFe))   #nutrient limitation for diaz
	
	dB_p =  (mu_p*gamma_p - m - kappa)*B_p                 #incremental change in non-diaz
	dB_d =  (mu_d*gamma_d - m - kappa)*B_d                 #incremental change in diaz
	dN   = -mu_p*gamma_p*B_p                                 + m*B_p         + m*B_d       + kappa*(N_0-N) #incremental change in N
	dP   = -mu_p*gamma_p*B_p*r_pP  - mu_d*gamma_d*B_d*r_dP*2 + m*B_p*r_pP*2  + m*B_d*r_dP  + kappa*(P_0-P)            #incremental change in P
	dFe  = -mu_p*gamma_p*B_p*r_pFe - mu_d*gamma_d*B_d*r_dFe  + m*B_p*r_pFe   + m*B_d*r_dFe + kappa*(Fe_0-Fe) + f_atm  #incremental change in Fe
    
    list(c(dB_p,dB_d,dN,dP,dFe))})}


times <- seq(0, 1000, by = 0.1)  #times for which we solve the system

out <- as.data.frame(ode(y=state, times=times, func=dxdt, parms=parameters))

out[nrow(out),2:6]

###########################################################################
## ANALYSIS ###############################################################
###########################################################################
N_f_atm <- 100
N_P0    <- 101
f_atms <- seq(0,0.0001,length.out=N_f_atm)
P_0s   <- seq(0.1*(1.6/16),2*(1.6/16),length.out=N_P0)

SS=SS_P=SS_N <- matrix(NA,ncol=N_f_atm,nrow=N_P0)
SS_pref=SS_P_pref=SS_N_pref <- matrix(NA,ncol=N_f_atm,nrow=N_P0)

for(i in 1:N_f_atm){
	print(i)
	parameters$f_atm <- f_atms[i]
	for(j in 1:N_P0){
		parameters$P_0 <- P_0s[j]
		
		out      <- as.data.frame(ode(y=state,times=times,func=dxdt,parms=parameters))
		out_pref <- as.data.frame(ode(y=state,times=times,func=dxdt_pref,parms=parameters)) 
		
		SS[j,i]   <- out[nrow(out),3]
		SS_P[j,i] <- out[nrow(out),5]
		SS_N[j,i] <- out[nrow(out),4]

		SS_pref[j,i]   <- out_pref[nrow(out),3]
		SS_P_pref[j,i] <- out_pref[nrow(out),5]
		SS_N_pref[j,i] <- out_pref[nrow(out),4]
	}
}

#######################################################################
## PLOTS ##############################################################
#######################################################################

par(mfrow=c(3,2),mar=c(4,4,4,5))


filled.contour3(y=(1/parameters$r_pP)*(P_0s/parameters$N_0),
			   x=(parameters$kappa*parameters$Fe_0 + f_atms)/(parameters$kappa*parameters$N_0)*(1/parameters$r_pFe),
			   z=t(SS_pref),
			   nlevels=40,
			   col=matlab.like(48),
			   zlim=c(0,1))
image.plot(matrix(seq(0,1,0.1)),legend.only=TRUE)

filled.contour3(y=(1/parameters$r_pP)*(P_0s/parameters$N_0),
			   x=(parameters$kappa*parameters$Fe_0 + f_atms)/(parameters$kappa*parameters$N_0)*(1/parameters$r_pFe),
			   z=t(SS_N),
			   nlevels=20,
			   col=matlab.like(20),
			   zlim=c(0,max(c(SS_N,SS_N_pref))))
filled.contour3(y=(1/parameters$r_pP)*(P_0s/parameters$N_0),
			   x=(parameters$kappa*parameters$Fe_0 + f_atms)/(parameters$kappa*parameters$N_0)*(1/parameters$r_pFe),
			   z=t(SS_N_pref),
			   nlevels=20,
			   col=matlab.like(20),
   			   zlim=c(0,max(c(SS_N,SS_N_pref))))
image.plot(matrix(c(0,max(c(SS_N,SS_N_pref)))),legend.only=TRUE)
			   
filled.contour3(y=(1/parameters$r_pP)*(P_0s/parameters$N_0),
			   x=(parameters$kappa*parameters$Fe_0 + f_atms)/(parameters$kappa*parameters$N_0)*(1/parameters$r_pFe),
			   z=t(SS_P),
			   nlevels=20,
			   col=matlab.like(20),
			   zlim=c(0,max(c(SS_P,SS_P_pref))))
filled.contour3(y=(1/parameters$r_pP)*(P_0s/parameters$N_0),
			   x=(parameters$kappa*parameters$Fe_0 + f_atms)/(parameters$kappa*parameters$N_0)*(1/parameters$r_pFe),
			   z=t(SS_P_pref),
			   nlevels=20,
			   col=matlab.like(20),
   			   zlim=c(0,max(c(SS_P,SS_P_pref))))
image.plot(matrix(c(0,max(c(SS_P,SS_P_pref)))),legend.only=TRUE)





