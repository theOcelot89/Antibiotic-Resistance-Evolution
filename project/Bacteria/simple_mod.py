import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def psi(a, psi_max, psi_min, zMIC, k):
    term = (a / zMIC)**k
    return psi_max - ((psi_max - psi_min) * term) / (term + 1)

def a(t):
    return 12.5 if t < 5 else 0

def dX_dt(X, t, psi_max, psi_min, zMIC, k):
    a_t = a(t)
    growth_rate = np.log(10) * psi(a_t, psi_max, psi_min, zMIC, k) * X
    return max(growth_rate, -X / 0.04)

# Model parameters
psi_min = -2
zMIC = 2
k = 0.8  # Using a single mean k value
psi_max = 0.8  # Assuming this is the psi_max for the simulation

# Time vector
t = np.linspace(0, 10, 101) ### changed 100 to 101 for even spacing

# Initial conditions
X0 = 1e8 ### tried with y as vector but didnt work

# Run a single simulation
X = odeint(dX_dt, X0, t, args=(psi_max, psi_min, zMIC, k))

# Plotting
fig, ax = plt.subplots()
ax.plot(t, X, label=f'Growth with k={k}')
ax.set_xlabel('Time')
ax.set_ylabel('Bacterial Density')
ax.set_yscale('log')
ax.set_ylim(1, 1e9)
ax.legend()
plt.savefig('simple_mod.png')
