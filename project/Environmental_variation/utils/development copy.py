import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the custom ODE function representing the variation
def variation(variation, t, dt, A, L, R):
    # Calculate the rate of change of variation
    d_variation_dt = A * np.sin(2 * np.pi * t / (L * R))
    
    # Update the variation by subtracting the previous value
    variation += d_variation_dt * dt
    
    return variation

# Define parameters
A = 1.0  # Amplitude
L = 1.0  # Wavelength
R = 1.0  # Frequency

# Time points
t = np.linspace(0, 10, 100)  # 100 time points from 0 to 10
dt = t[1] - t[0]  # Time step size

# Initial condition
variation0 = 0.0  # Initial value of the variation

# Solve the ODE
variation_over_time = odeint(variation, variation0, t, args=(dt, A, L, R))

# Plot the variation
plt.plot(t, variation_over_time)
plt.xlabel('Time')
plt.ylabel('Variation')
plt.title('Variation Over Time')
plt.grid(True)
plt.show()
