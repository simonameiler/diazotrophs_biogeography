import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# Euler forward
n  = 10000
N = np.zeros(n) 
P = np.zeros(n)

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

###################################################
## NUMERICAL SOLUTION #############################
###################################################
# specify indices
i_p = 0
i_n = 1

# specify parameter values
theta = {
    'mu':1.0,
    'K':0.1,
    'kappa':0.05,
    'm':0.1,
    'n0':0.2,
}

# specify model 
def dxdt(x,t,theta):
    growth = theta['mu']*(x[i_n]/(x[i_n] + theta['K']))*x[i_p]
    p_mort = theta['m']*x[i_p]
    p_mix  = theta['kappa']*x[i_p]
    n_mix  = theta['kappa']*(theta['n0'] - x[i_n])
    dP     = growth - p_mort - p_mix
    dN     = p_mort -growth + n_mix
    return np.array((dP, dN))

# initial conditions
x0 = (1.0,1.0)

# times where you want the solution
t = np.arange(0.0,1000.0,0.01)

# solve ODE
x = odeint(dxdt, x0, t, args=(theta,))

# plot
fig,ax = plt.subplots()
colors = ('C1','C2')

for i,name in enumerate(('P','N')):
    ax.plot(t,x[:,i], color=colors[i], label=name)

ax.legend()

