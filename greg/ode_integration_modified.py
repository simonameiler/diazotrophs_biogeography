import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import cmocean as cm


# Euler forward
n  = 10000
N = np.zeros(n) #nurtient
P = np.zeros(n) #

mu    = 0.01
K     = 0.1
kappa = 0.05
m     = 0.1
n_0   = 0.05

N[0] = 1
P[0] = 1

dt = 0.1



for t in range(1,n):
    growth = mu*(N[t-1]/(N[t-1]+K))*P[t-1]
    p_mort = m*P[t-1]
    p_mix  = kappa*P[t-1]
    n_mix  = kappa*(n_0 - N[t-1])
    
    dP = (growth - p_mort - p_mix)*dt
    dN = (p_mort - growth + n_mix)*dt
    
    N[t] = N[t-1] + dN
    P[t] = P[t-1] + dP
    
fig,ax = plt.subplots()
colors = ('C1','C2')

ax.plot(N, color='C1', label='N')
ax.plot(P, color=colors[1], label='P')
ax.legend()

#%%################################################
## NUMERICAL SOLUTION #############################
###################################################
# specify indices
i_p  = 0
i_d  = 1  
i_N  = 2
i_P  = 3 
i_Fe = 4

# specify parameter values
theta = {
    'mu_p':2.5,         #d-1
    'mu_d':1.25,        #d-1
    'K_p_N':0.056,      #mmol N m-3
    'K_p_P':0.035,      #mmol P m-3
    'K_p_Fe':0.00035,   #mmol Fe m-3
    'K_d_P':0.035,      #mmol P m-3
    'K_d_Fe':0.0011,    #mmol Fe m-3
    'r_p_P':0.0625,     # -
    'r_p_Fe':6.25e-05,  # -
    'r_d_P':0.025,      # -
    'r_d_Fe':7.5e-04,   # -
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
    d_growth = theta['mu_d']*min((x[i_P]/(x[i_P]+theta['K_d_P'])),(x[i_Fe]/(x[i_Fe]+theta['K_d_Fe'])))*x[i_d]
    p_mort = theta['m']*x[i_p]
    d_mort = theta['m']*x[i_d]
    p_mix  = theta['kappa']*x[i_p]
    d_mix  = theta['kappa']*x[i_d]
    N_mix  = theta['kappa']*(theta['N0'] - x[i_N])
    P_mix  = theta['kappa']*(theta['P0'] - x[i_P])
    Fe_mix = theta['kappa']*(theta['Fe0'] - x[i_Fe])
    dp     = p_growth - p_mort - p_mix
    dd     = d_growth - d_mort - d_mix
    dN     = -p_growth + p_mort + d_mort + N_mix
    dP     = -p_growth*theta['r_p_P'] - d_growth*theta['r_d_P'] + p_mort*theta['r_p_P'] + d_mort*theta['r_d_P'] + P_mix
    dFe    = -p_growth*theta['r_p_Fe'] - d_growth*theta['r_d_Fe'] + p_mort*theta['r_p_Fe'] + d_mort*theta['r_d_Fe'] + Fe_mix + theta['f_atm']
    return np.array((dp, dd, dN, dP, dFe))

# initial conditions
x0 = (1.0,1.0,1.0,1.0,1.0)

# times where you want the solution
t = np.arange(0.0,1000.0,0.01)

# solve ODE
x = odeint(dxdt, x0, t, args=(theta,))

# plot
fig,ax = plt.subplots()
colors = ('C0','C1','C2','C3','C4')

for i,name in enumerate(('p','d','N','P','Fe')):
    ax.plot(t,x[:,i], color=colors[i], label=name)

ax.legend()

#%% Modify the simple framework of above
# initial conditions
x0 = (1.0,1.0,1.0,1.0,1.0)

# times where you want the solution
t = np.arange(0.0,1000.0,0.01)

# build loop here to vary P0 and f_atm
var_P0 = np.linspace(0.08,0.16,50)
var_f_atm = np.linspace(1e-05,1e-01,50)
phi_PN = np.zeros(len(var_P0))
phi_FeN = np.zeros(len(var_f_atm))
matrix_steady_state = np.zeros((len(var_P0),len(var_f_atm),5))
for i in range(len(var_P0)):
    theta['P0'] = var_P0[i]
    phi_PN[i] = theta['P0']/theta['N0']*(1/theta['r_p_P'])
    for j in range(len(var_f_atm)):
        theta['f_atm'] = var_f_atm[j]
        phi_FeN[j] = ((theta['kappa']*theta['Fe0']+theta['f_atm'])/theta['N0']*theta['kappa'])*(1/theta['r_p_Fe'])
        x = odeint(dxdt, x0, t, args=(theta,))
        #print('i: '+str(i))
        #print('j: '+str(j))
        #print('f_atm: '+str(theta['f_atm']))
        #print('P0: '+str(theta['P0']))
        #print(x[-1,:])
        matrix_steady_state[i,j,:] = x[-1,:]


# solve ODE
#x = odeint(dxdt, x0, t, args=(theta,))

#%% plot
fig,ax = plt.subplots()
colors = ('C0','C1','C2','C3','C4')

for i,name in enumerate(('p','d','N','P','Fe')):
    ax.plot(t,x[:,i], color=colors[i], label=name)

ax.legend()

#%% 
fig,ax = plt.subplots(1,2,figsize=(9,4))#,sharey=True)
c0 = ax[0].contourf(phi_FeN, phi_PN, matrix_steady_state[:,:,4],cmap=cm.cm.haline,levels=np.linspace(0,1,101),extend='both')
c1 = ax[1].contourf(phi_FeN, phi_PN, matrix_steady_state[:,:,1],cmap=cm.cm.haline,levels=np.linspace(0,2,101),extend='both')
for i in range(0,2):
    #ax[i].axhline(var_P0,linewidth=1.0,linestyle='dashed',color='w')
    #ax[i].axvline(var_f_atm,linewidth=1.0,linestyle='dashed',color='w')
    ax[i].set_ylabel('P:N')
    ax[i].set_xlabel('Fe:N')
#ax[0].text(0.9,0.95,'area',transform=ax[0].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#ax[1].text(0.85,0.95,'accuracy',transform=ax[1].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
cbar0 = plt.colorbar(c0,ax=ax[0])
#cbar0.set_label('(m)',rotation=90, position=(0.5,0.5))
cbar1 = plt.colorbar(c1,ax=ax[1])
#cbar1.set_label('(-)',rotation=90, position=(0.5,0.5))
plt.show()


#%% Modify and visualize other parameters

# build loop here to vary P0 and f_atm
x0 = (0.1,0.1,0.1,0.1,0.1)

var_rdFe = [3.75e-04, 7.5e-04, 1.5e-03]
var_f_atm = np.logspace(-6,-3,50)

theta['P0'] = 0.12
phi_PN = theta['P0']/theta['N0']*(1/theta['r_p_P'])

phi_FeN = np.zeros(len(var_f_atm))
matrix_steady_state = np.zeros((len(var_rdFe),len(var_f_atm),5))
for i in range(len(var_rdFe)):
    theta['r_d_Fe'] = var_rdFe[i]
    for j in range(len(var_f_atm)):
        theta['f_atm'] = var_f_atm[j]
        phi_FeN[j] = ((theta['kappa']*theta['Fe0']+theta['f_atm'])/theta['N0']*theta['kappa'])*(1/theta['r_p_Fe'])
        x = odeint(dxdt, x0, t, args=(theta,))
        #print('i: '+str(i))
        #print('j: '+str(j))
        #print('f_atm: '+str(theta['f_atm']))
        #print('P0: '+str(theta['P0']))
        #print(x[-1,:])
        matrix_steady_state[i,j,:] = x[-1,:]
        
        
#%% plot diazotroph vs. phi_FeN
fig,ax = plt.subplots(1,1,figsize=(9,4))#,sharey=True)
plt.plot(phi_FeN[:],matrix_steady_state[0,:,0],label='rdFe= 0.5x')
plt.plot(phi_FeN[:],matrix_steady_state[1,:,0],label='rdFe= 1x')
plt.plot(phi_FeN[:],matrix_steady_state[2,:,0],label='rdFe= 2x')
plt.plot(phi_FeN[:],matrix_steady_state[0,:,1],label='rdFe= 0.5x')
plt.plot(phi_FeN[:],matrix_steady_state[1,:,1],label='rdFe= 1x')
plt.plot(phi_FeN[:],matrix_steady_state[2,:,1],label='rdFe= 2x')
ax.legend(loc='upper right')
ax.set_xlabel('normalized Fe:N supply ratio')
ax.set_ylabel('diazotrophs (mmol N m-3)')
plt.show()