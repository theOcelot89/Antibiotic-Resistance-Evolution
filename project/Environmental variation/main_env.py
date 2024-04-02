import numpy as np
import matplotlib.pyplot as plt

# Define the function to model environmental variation
def environmental_variation(t, A, B, L, R, epsilon):
    return A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon

# Set common parameters
T = 1000 # number of time steps to simulate (e.g., days in a year)
L = [5] # list for different lifespans
epsilon = np.random.normal(0, 1, T) # stochastic error term

# Parameters to vary
A_values = [0, 1] # Amplitude of deterministic component
B_values = [0.5, 0.1] # Amplitude of stochastic component
R_value = 10 # scale of 

for L in L:
    # Create plots
    fig, axs = plt.subplots(2, 2, figsize=(15, 20)) #creates an empty figure with 2x2 grid
    t = np.arange(T)

    for i, A in enumerate(A_values):
        # Varying A, keeping B, L, and R constant
        E = environmental_variation(t, A, 1, L, R_value, epsilon)
        axs[i, 0].plot(t, E)
        axs[i, 0].set_title(f'Varying A (A={A}, B=1, R={R_value}), lifespan = {L}')
        axs[i, 0].set_xlabel('Time')
        axs[i, 0].set_ylabel('Environmental Variation')

    for i, B in enumerate(B_values):
        # Varying B, keeping A, L, and R constant
        E = environmental_variation(t, 1, B, L, R_value, epsilon)
        axs[i, 1].plot(t, E)
        axs[i, 1].set_title(f'Varying B (A=1, B={B}, R={R_value}), lifespan = {L}')
        axs[i, 1].set_xlabel('Time')
        axs[i, 1].set_ylabel('Environmental Variation')


    plt.tight_layout()
    plt.savefig(f'lifespan = {L} environmental_variation_gridh.jpeg')








