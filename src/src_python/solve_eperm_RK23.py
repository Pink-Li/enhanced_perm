import numpy as np
from scipy.integrate import solve_ivp
#import matplotlib.pyplot as plt
import timeit

nx = 250  # number of grid points
L = 2.5  #length of domain
X = np.linspace(0, L, nx) # position along the domain
h = L / nx #grid spacing
t=25 #time period of solving
nt=250 #time steps
alpha = 0.001 #diffusivity
gamma=0.05 #pressure sensitivity
u_0=100.0 #BC1
u_end=1.0 #BC2

def odefunc(t, u):
    dudt = np.zeros(X.shape)

    dudt[0] = 0 # no flux boundary condition
    dudt[-1] = 0

    # now for the internal nodes
    for i in range(1, nx-1):
        dudt[i] = alpha * gamma*u[i]* (u[i + 1] - 2*u[i] + u[i - 1]) / h**2
    return dudt

init = np.exp(u_end * np.ones(X.shape)*gamma)/gamma # initializing the pressure vector
init[0] = np.exp(u_0*gamma)/gamma  # cole-hopf transformation at one boundary condition
init[-1] = np.exp(u_end*gamma)/gamma #cole-hopf transformation the other boundary condition
tspan=(0,t)
teval = np.linspace(0, t, nt)
start=timeit.default_timer()
sol = solve_ivp(odefunc, tspan, init, method='RK23',t_eval=teval)
stop=timeit.default_timer()
sol_use=np.transpose(sol.y)
pxt=np.log(sol_use*gamma)/gamma
comp_time=stop-start
