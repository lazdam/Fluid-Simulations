"""
Evolving an adiabatic shock wave
@author: Mattias Lazda
@collab: Jules Fauscher 
March 19th 2022
"""
import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing
Ngrid = 100
Nsteps = 5000
dt = 0.1
dx = 2.0
gamma = 5/3 # adiabatic index

x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
f3 = np.ones(Ngrid) # rho x e_tot
u = np.zeros(Ngrid+1) # advective velocity (keep the 1st and last element zero)


# Pressure 
P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)

# Compute sound speed squared
cs2 = gamma * P / f1

# Get initial mach number
u_wave = f2/f1
mach = u_wave/np.sqrt(cs2)

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update

    return f

# Apply initial Gaussian perturbation to energy 
Amp, sigma = 0.5, Ngrid/10
f2 = f2 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
f1 = f1 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
f3 = f3 + 1000 * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

# plotting
pl.ion()
fig, ax = pl.subplots(2,1)

# Density plot
x1, = ax[0].plot(x, f1, 'ro')
#ax[0].set_xlim([0, dx*Ngrid+1])
#ax[0].set_ylim([0, 10])
ax[0].set_xlabel('x')
ax[0].set_ylabel('Density')

# Mach number plot
x2, = ax[1].plot(x, mach, 'ro')
#ax[1].set_xlim([0, dx*Ngrid+1])
ax[1].set_ylim([-10, 10])
ax[1].set_xlabel('x')
ax[1].set_ylabel(r'$\mathcal{M}$')

fig.canvas.draw()

for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:])) # Need to add third velocity 

    # update density and momentum
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)

    # Compute pressure 
    P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)

    # Get cs2 
    cs2 = gamma*P/f1
    
    # Add pressure gradient 
    f2[1:-1] = f2[1:-1] - (dt/dx) * cs2[1:-1] * (f1[2:] - f1[:-2])

    # Correct source term at boundary 
    f2[0] = f2[0] - 0.5 * (dt/dx) * cs2[0] * (f1[1] - f1[0])
    f2[-1] = f2[-1] - 0.5 * (dt/dx) * cs2[-1] * (f1[-1] - f1[-2])

    
    # Re-calculate advection velocities
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:])) 

    # Advect energy 
    f3 = advection(f3, u, dt, dx)

    # Re calculate pressure
    P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)

    # Add source to energy term
    f3[1:-1] = (gamma/(gamma-1))*P[1:-1] + 0.5*(f2[1:-1]**2/f1[1:-1]) 

    #correct for source term at the boundary (reflective)
    f3[0] = (gamma/(gamma-1))*P[0] + 0.5*(f2[0]**2/f1[0])
    f3[-1] = (gamma/(gamma-1))*P[0] + 0.5*(f2[0]**2/f1[0])

    # Update pressure and sound speed
    P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)
    cs2 = gamma * P / f1
    
    u_wave = f2/f1
    mach = u_wave/np.sqrt(cs2)

    # update the plot
    x1.set_ydata(f1)
    x2.set_ydata(mach)
    fig.canvas.draw()
    pl.pause(0.01)

