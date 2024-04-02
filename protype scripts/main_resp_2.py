import numpy as np
import matplotlib.pyplot as plt

# Define the function to model environmental variation and reaction norm
def environmental_variation(t, A, B, L, R, epsilon):
    return A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon

def reaction_norm(I0, b, C):
    return I0 + b * C

# Set common parameters for environmental variation, adjusting B to introduce stochasticity
T = 365 # number of time steps to simulate (e.g., days in a year)
L = 5 # lifespan
epsilon = np.random.normal(0, 1, T) # stochastic error term
A, B, L, R = 1, 0.2, 365, 0.1 # Adjusted B to include stochastic error in environmental variation
E = environmental_variation(np.arange(T), A, B, L, R, epsilon)

# Parameters for three genotypes
genotypes = {
    "Genotype 1": {"I0": 0, "b": 0.5},
    "Genotype 2": {"I0": 0, "b": 1},
    "Genotype 3": {"I0": 0.5, "b": 0.04}
}
colors = ['blue', 'green', 'red'] # Define colors for each genotype

# Plot the reaction norms for genotypes across the environmental conditions
plt.figure(figsize=(10, 6))

# Plot environmental variation with stochasticity
plt.plot(np.arange(T), E, label='Environmental Variation', linestyle='--', color='grey')

plt.savefig('environmental_variation.png')

for (name, params), color in zip(genotypes.items(), colors):
    I = reaction_norm(params["I0"], params["b"], E)
    plt.plot(np.arange(T), I, label=name, color=color)

plt.title('Reaction Norms Across Time')
plt.xlabel('Time')
plt.ylabel('Phenotypic Expression (I)')
plt.legend()

# Save the plot to a file
plt.savefig('Reaction Norms Across Time.png')

# Prepare for the second plot
plt.figure(figsize=(10, 6))
colors = ['blue', 'green', 'red'] # Define colors for each genotype

for (name, params), color in zip(genotypes.items(), colors):
    I = reaction_norm(params["I0"], params["b"], E)
    plt.plot(E, I, label=name, color=color)

plt.title('Reaction Norms for Genotypes')
plt.xlabel('Environmental Variation (E)')
plt.ylabel('Phenotypic Expression (I)')
plt.legend()

plt.savefig('Reaction Norms for Genotypes.png') # Save the plot