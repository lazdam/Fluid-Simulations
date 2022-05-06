import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.family':'serif'})
# PART 1.2
def cost(u, m):

    U = (u/10)
    M = (m/0.035)
    term1 = 10.7/U * M**(-1/4)
    term2 = 2.3e-4 * U**2 * M**(-1/3) 
    term3 = 2.5/U**2 * M**(1/3)

    return term1 + term2 + term3
    
# Define velocities
u = np.linspace(0.1, 1000, 10000)

# Get log spaced values
mass = np.logspace(-2, 2, 5)

# Store minimum values of u for each mass
min_u = []
min_eps = []
plt.figure(figsize=(15,5))
for i in mass: 
    eps = cost(u, i)

    # Save minimum u and minimum eps for later
    min_u.append(eps.argmin())
    min_eps.append(np.min(eps))

    # Plot
    plt.plot(u, eps, label = 'm = {0} kg'.format(i))

plt.xlabel(r'Velocity, $u$ (m/s)')
plt.ylabel(r'$\epsilon$')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.grid()
plt.title(r'$\epsilon$ as a function of $u$')
plt.savefig('Q1_2.png')
plt.show()


# PART 1.3

# Determine location of minimum value
min_u = np.asarray(min_u)
min_eps = np.asarray(min_eps)

def power_law(u, a, b):
    return a*u**b 

xx = np.logspace(-2,2,1000)
p0 = np.array([2,-0.4])
yy = power_law(xx,*p0)

fit, cov = curve_fit(power_law, mass, min_eps, p0=p0)
errs = np.sqrt(np.diag(cov))
y_fit = power_law(xx, *fit)

plt.figure(figsize=(15,5))
plt.plot(mass, min_eps, 'o', label = 'Data', color = 'black')
plt.plot(xx, y_fit, label = 'Fit', color = 'firebrick')
plt.legend()
plt.xscale('log')
plt.xlabel('Mass, m (kg)')
plt.ylabel(r'$\epsilon_{min}$')
plt.grid()
plt.text(10, 0.3, r'$\epsilon_{min} = a\times m^b$', fontsize='14')
plt.text(10, 0.25, r'$a = {0} \pm {1}$'.format(round(fit[0],3), round(errs[0],3)), fontsize='12')
plt.text(10, 0.20, r'$b = {0} \pm {1}$'.format(round(fit[1],3), round(errs[1],3)), fontsize='12')
plt.title(r'$\epsilon$ as a function of $m$')
plt.savefig('Q1_3.png')
plt.show()

