"""
Evolving an adiabatic shock wave
@author: Mattias Lazda
@collab: Jules Fauscher, Nicolas Desjardins, David Tucci
March 19th 2022
"""
import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing
Ngrid = 100
Nsteps = 5000
dt = 0.012
dx = 1.5
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
Amp, sigma = 10000, Ngrid/10
f3 = f3 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

# Get analytic solution
rho2 = np.empty(Ngrid)
rho2.fill((gamma + 1)/(gamma-1))

# plotting
pl.ion()
fig, ax = pl.subplots(2,1)

# Density plot
x1, = ax[0].plot(x, f1, 'ro', color='firebrick')
x11, = ax[0].plot(x, rho2, label = 'Strong Shock Solution', color='blue')
ax[0].set_xlim([0, dx*Ngrid])
ax[0].set_xlim([0, dx*Ngrid])
ax[0].set_ylim([0, 5])
ax[0].grid()
ax[0].set_xlabel('x')
ax[0].set_ylabel('Density')
ax[0].legend(loc='upper center')

# Mach number plot
x2, = ax[1].plot(x, mach, 'ro', color='firebrick')
ax[1].set_xlim([0, dx*Ngrid])
ax[1].set_ylim([-2, 2])
ax[1].set_xlabel('x')
ax[1].grid()
ax[1].set_ylabel(r'$\mathcal{M}$')

fig.canvas.draw()

for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:])) 

    # update density and momentum
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)

    # Compute pressure 
    P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)

    # Get cs2 
    cs2 = gamma*P/f1
    
    # Add pressure gradient 
    f2[1:-1] = f2[1:-1] -  0.5*(dt/dx) * (P[2:] - P[:-2]) 

    # Correct source term at boundary 
    f2[0] = f2[0] -  0.5*(dt/dx) * (P[1] - P[0]) 
    f2[-1] = f2[-1] -  0.5*(dt/dx) * (P[-1] - P[-2])

    # Re-calculate advection velocities
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:])) 

    # Advect energy 
    f3 = advection(f3, u, dt, dx)

    # Re calculate pressure
    P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)

    # The below term is needed to evolve energy in time
    u_wave = f2/f1
    Pu = P*u_wave

    # Add source to energy term
    f3[1:-1] = f3[1:-1] - 0.5*(dt/dx)*(Pu[2:] - Pu[:-2]) 

    #correct for source term at the boundary (reflective)
    f3[0] = f3[0] -  0.5*(dt/dx)*(Pu[1] - Pu[0]) 
    f3[-1] = f3[-1] -  0.5*(dt/dx)*(Pu[-1]-Pu[-2])

    # Update pressure and sound speed
    P = (gamma - 1)/gamma * (f3 - 0.5*f2**2/f1)
    cs2 = gamma * P / f1
    
    # Calculate wave speed and mach number
    u_wave = f2/f1
    mach = u_wave/np.sqrt(cs2)



    # Update the plot
    x1.set_ydata(f1)
    x2.set_ydata(mach)
    fig.canvas.draw()
    pl.pause(0.01)

