"""
PHYS432 W2022
Diffusion
@author: Mattias Lazda
@collab: Nicholas Desjardins, Jules Fauscher, David Tucci
March 9, 2022

NOTE: This code is heavily based off of the sample code provided by Prof. Lee, available here: https://github.com/evjhlee/phys432-w2022/blob/main/phys432_w2022_advdiff.py
"""
import numpy as np
import matplotlib.pyplot as pl

# Set up the grid and advection and diffusion parameters
Ngrid = 120
Nsteps = 10000
dt = 9
dx = 1
g = 0.01

# Set the diffusion coefficient
D = 0.5

# Set arbitrary tilt to the plane
alpha = 10 #deg 
a = g*np.sin(np.radians(alpha))

# Courant condition is satisfied
beta = D*dt/dx**2

x = np.arange(0, Ngrid, dx, dtype = 'float') 

f = np.zeros(len(x))

# Set up plot
pl.ion()
fig, ax = pl.subplots(1,1)

# Plot analytic solution
u_analytic = -a/D * (0.5*np.power(x,2) - Ngrid*x)
ax.plot(x, u_analytic, label = 'Analytic Solution', color = 'darkorange')


# We will be updating these plotting objects
plt, = ax.plot(x, f, color = 'firebrick', label = 'Numerical Solution', ls = '--')
ax.legend()
ax.set_xlabel('Vertical position, arbitrary')
ax.set_ylabel('Velocity, arbitrary')

# Draw
fig.canvas.draw()

npts = len(x)
for ct in range(Nsteps):

    ## Calculate diffusion first
    # Setting up matrices for diffusion operator
    A = np.eye(npts) * (1.0 + 2.0 * beta) + np.eye(npts, k=1) * -beta + np.eye(npts, k=-1) * -beta

    # Apply no-slip on left edge 
    A[0] = np.zeros(npts)
    A[0][0] = 1.0
    
    # Stress-free boundary condition on the right
    A[-1][-1] = 1.0 + beta

    # Solving for the next timestep
    f = np.linalg.solve(A, f)
    # Take into account gravity
    f += a*dt    
    f[0] = 0

    # only plot every 10th step to speed up plot
    if ct%10==0: 
        # update the plot
        plt.set_ydata(f)

        fig.canvas.draw()
        ax.set_title('Step {0}'.format(ct))
        pl.pause(0.01)