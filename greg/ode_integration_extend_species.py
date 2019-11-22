import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import cmocean as cm

#%%################################################
## NUMERICAL SOLUTION #############################
###################################################
# specify indices
i_p  = 0
i_d1 = 1
i_d2 = 2  
i_d3 = 3
i_d4 = 4   
i_N  = 5
i_P  = 6 
i_Fe = 7

# specify parameter values - 4 different diazotroph species
theta = {
    'mu_p':2.5,         #d-1
    'mu_d1':1.25,       #d-1
    'mu_d2':1.0,        #d-1
    'mu_d3':1.5,        #d-1
    'mu_d4':1.75,       #d-1
    'K_p_N':0.056,      #mmol N m-3
    'K_p_P':0.035,      #mmol P m-3
    'K_p_Fe':0.00035,   #mmol Fe m-3
    'K_d1_P':0.035,     #mmol P m-3
    'K_d1_Fe':0.0011,   #mmol Fe m-3
    'K_d2_P':0.035,     #mmol P m-3
    'K_d2_Fe':0.0011,   #mmol Fe m-3
    'K_d3_P':0.035,     #mmol P m-3
    'K_d3_Fe':0.0011,   #mmol Fe m-3
    'K_d4_P':0.035,     #mmol P m-3
    'K_d4_Fe':0.0011,   #mmol Fe m-3
    'r_p_P':0.0625,     # -
    'r_p_Fe':6.25e-05,  # -
    'r_d1_P':0.025,     # -
    'r_d1_Fe':7.5e-04,  # -
    'r_d2_P':0.05,      # -
    'r_d2_Fe':7.5e-04,  # -
    'r_d3_P':0.001,     # -
    'r_d3_Fe':7.5e-03,  # -
    'r_d4_P':0.1,       # -
    'r_d4_Fe':7.5e-03,  # -
    'kappa':0.1,        #d-1
    'm':0.05,           #d-1
    'N0':1.6,           #mmol N m-3
    'P0':0.2,           #mmol P m-3 - variable?
    'Fe0':1e-05,        #mmol Fe m-3
    'f_atm':1e-03,      #mmol Fe m-3 d-1 - assumption from Ward et al. Fig. 6
}

# specify model 
def dxdt(x,t,theta):
    p_growth = theta['mu_p']*min((x[i_N]/(x[i_N]+theta['K_p_N'])),(x[i_P]/(x[i_P]+theta['K_p_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_p_Fe'])))*x[i_p]
    d1_growth = theta['mu_d1']*min((x[i_P]/(x[i_P]+theta['K_d1_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d1_Fe'])))*x[i_d1]
    d2_growth = theta['mu_d2']*min((x[i_P]/(x[i_P]+theta['K_d2_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d2_Fe'])))*x[i_d2]
    d3_growth = theta['mu_d3']*min((x[i_P]/(x[i_P]+theta['K_d3_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d3_Fe'])))*x[i_d3]
    d4_growth = theta['mu_d4']*min((x[i_P]/(x[i_P]+theta['K_d4_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d4_Fe'])))*x[i_d4]
    p_mort = theta['m']*x[i_p]
    d1_mort = theta['m']*x[i_d1]
    d2_mort = theta['m']*x[i_d2]
    d3_mort = theta['m']*x[i_d3]
    d4_mort = theta['m']*x[i_d4]
    p_mix  = theta['kappa']*x[i_p]
    d1_mix  = theta['kappa']*x[i_d1]
    d2_mix  = theta['kappa']*x[i_d2]
    d3_mix  = theta['kappa']*x[i_d3]
    d4_mix  = theta['kappa']*x[i_d4]
    N_mix  = theta['kappa']*(theta['N0'] - x[i_N])
    P_mix  = theta['kappa']*(theta['P0'] - x[i_P])
    Fe_mix = theta['kappa']*(theta['Fe0'] - x[i_Fe])
    dp     = p_growth - p_mort - p_mix
    dd1    = d1_growth - d1_mort - d1_mix
    dd2    = d2_growth - d2_mort - d2_mix
    dd3    = d3_growth - d3_mort - d3_mix
    dd4    = d4_growth - d4_mort - d4_mix
    dN     = -p_growth + p_mort + N_mix
    dP     = -p_growth*theta['r_p_P'] - d1_growth*theta['r_d1_P'] - d2_growth*theta['r_d2_P'] - d3_growth*theta['r_d3_P'] - d4_growth*theta['r_d4_P'] + theta['m']*(x[i_p]*theta['r_p_P']) + theta['m']*(x[i_d1]*theta['r_d1_P']) + theta['m']*(x[i_d2]*theta['r_d2_P']) + theta['m']*(x[i_d3]*theta['r_d3_P']) + theta['m']*(x[i_d4]*theta['r_d4_P']) + P_mix
    dFe    = -p_growth*theta['r_p_Fe'] - d1_growth*theta['r_d1_Fe'] - d2_growth*theta['r_d2_Fe'] - d3_growth*theta['r_d3_Fe'] - d4_growth*theta['r_d4_Fe'] + theta['m']*(x[i_p]*theta['r_p_Fe']) + theta['m']*(x[i_d1]*theta['r_d1_Fe']) + theta['m']*(x[i_d2]*theta['r_d2_Fe']) + theta['m']*(x[i_d3]*theta['r_d3_Fe']) + theta['m']*(x[i_d4]*theta['r_d4_Fe']) + Fe_mix + theta['f_atm']
    return np.array((dp, dd1, dd2, dd3, dd4, dN, dP, dFe))

# initial conditions
x0 = (1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)

# times where you want the solution
t = np.arange(0.0,1000.0,0.01)

# solve ODE
x = odeint(dxdt, x0, t, args=(theta,))

# plot
fig,ax = plt.subplots()
colors = ('C0','C1','C2','C3','C4','C5','C6','C7')

for i,name in enumerate(('p','d1','d2','d3','d4','N','P','Fe')):
    ax.plot(t,x[:,i], color=colors[i], label=name)

ax.legend()


#%% specify indices
i_p  = 0
i_d1 = 1
i_d2 = 2    
i_N  = 3
i_P  = 4 
i_Fe = 5

# specify parameter values - 2 different diazotroph species
theta = {
    'mu_p':2.5,         #d-1
    'mu_d1':1.25,       #d-1
    'mu_d2':1.25,        #d-1
    'K_p_N':0.056,      #mmol N m-3
    'K_p_P':0.035,      #mmol P m-3
    'K_p_Fe':0.00035,   #mmol Fe m-3
    'K_d1_P':0.035,     #mmol P m-3
    'K_d1_Fe':0.0011,   #mmol Fe m-3
    'K_d2_P':0.035,     #mmol P m-3
    'K_d2_Fe':0.0011,   #mmol Fe m-3
    'r_p_P':0.0625,     # -
    'r_p_Fe':6.25e-05,  # -
    'r_d1_P':0.025,     # -
    'r_d1_Fe':7.5e-04,  # -
    'r_d2_P':0.025,      # -
    'r_d2_Fe':7.5e-04,  # -
    'kappa':0.1,        #d-1
    'm':0.05,           #d-1
    'N0':1.6,           #mmol N m-3
    'P0':0.2,           #mmol P m-3 - variable?
    'Fe0':1e-05,        #mmol Fe m-3
    'f_atm':1e-03,      #mmol Fe m-3 d-1 - assumption from Ward et al. Fig. 6
}

# specify model 
def dxdt(x,t,theta):
    p_growth = theta['mu_p']*min((x[i_N]/(x[i_N]+theta['K_p_N'])),(x[i_P]/(x[i_P]+theta['K_p_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_p_Fe'])))*x[i_p]
    d1_growth = theta['mu_d1']*min((x[i_P]/(x[i_P]+theta['K_d1_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d1_Fe'])))*x[i_d1]
    d2_growth = theta['mu_d2']*min((x[i_P]/(x[i_P]+theta['K_d2_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d2_Fe'])))*x[i_d2]
    p_mort = theta['m']*x[i_p]
    d1_mort = theta['m']*x[i_d1]
    d2_mort = theta['m']*x[i_d2]
    p_mix  = theta['kappa']*x[i_p]
    d1_mix  = theta['kappa']*x[i_d1]
    d2_mix  = theta['kappa']*x[i_d2]
    N_mix  = theta['kappa']*(theta['N0'] - x[i_N])
    P_mix  = theta['kappa']*(theta['P0'] - x[i_P])
    Fe_mix = theta['kappa']*(theta['Fe0'] - x[i_Fe])
    dp     = p_growth - p_mort - p_mix
    dd1    = d1_growth - d1_mort - d1_mix
    dd2    = d2_growth - d2_mort - d2_mix
    dN     = -p_growth + p_mort + N_mix
    dP     = -p_growth*theta['r_p_P'] - d1_growth*theta['r_d1_P'] - d2_growth*theta['r_d2_P'] + theta['m']*(x[i_p]*theta['r_p_P']) + theta['m']*(x[i_d1]*theta['r_d1_P']) + theta['m']*(x[i_d2]*theta['r_d2_P']) + P_mix
    dFe    = -p_growth*theta['r_p_Fe'] - d1_growth*theta['r_d1_Fe'] - d2_growth*theta['r_d2_Fe'] + theta['m']*(x[i_p]*theta['r_p_Fe']) + theta['m']*(x[i_d1]*theta['r_d1_Fe']) + theta['m']*(x[i_d2]*theta['r_d2_Fe']) + Fe_mix + theta['f_atm']
    return np.array((dp, dd1, dd2, dN, dP, dFe))

# initial conditions
x0 = (1.0,1.0,1.0,1.0,1.0,1.0)

# times where you want the solution
t = np.arange(0.0,1000.0,0.01)

# solve ODE
x = odeint(dxdt, x0, t, args=(theta,))

# plot
fig,ax = plt.subplots()
colors = ('C0','C1','C2','C3','C4','C5')

for i,name in enumerate(('p','d1','d2','N','P','Fe')):
    ax.plot(t,x[:,i], color=colors[i], label=name)

ax.legend()

