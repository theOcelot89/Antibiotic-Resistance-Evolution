import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the differential equations with mutation
def population_dynamics(y, t, r1, r2, t_mutation):
    N = y[0]
    
    if t < t_mutation:
        dNdt = r1 * N  # Before mutation
    else:
        dNdt = r2 * N  # After mutation
    
    return [dNdt]

# Initial condition
N0 = [1]  # Initial bacteria population

# Time points where solution is computed
t = np.linspace(0, 10, 100)

# Parameters
r1 = 0.1       # Growth rate before mutation
r2 = 0.2       # Growth rate after mutation
t_mutation = 5 # Time of mutation

# Solve the ODE
solution = odeint(population_dynamics, N0, t, args=(r1, r2, t_mutation))

# Plot the results
plt.plot(t, solution[:, 0], label='Bacteria Population')
plt.axvline(x=t_mutation, color='r', linestyle='--', label='Mutation Time')
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.show()
