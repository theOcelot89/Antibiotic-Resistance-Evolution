import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the function to model environmental variation and reaction norm
def environmental_variation(t, A, B, L, R, epsilon):
    return A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon

def reaction_norm(I0, b, C):
    return I0 + b * C

# Set common parameters for environmental variation, adjusting B to introduce stochasticity
T = 110 # number of time steps to simulate (e.g., days in a year)
epsilon = np.random.normal(0, 1, T) # stochastic error term
A, B, L, R = 1, 0.1, 10, 100 # Adjusted B to include stochastic error in environmental variation
E = environmental_variation(np.arange(T), A, B, L, R, epsilon)


# if any of the values are negative, set them to 0
E[E < 0] = 0.1

# if any of the values are greater than 1, set them to 1
E[E > 1] = 1

# Parameters for three genotypes
genotypes = {
    "Genotype 1": {"I0": 0, "b": 0.5},
    "Genotype 2": {"I0": 0, "b": 1},
    "Genotype 3": {"I0": 0.5, "b": 0}
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


def psi(a, psi_max, psi_min, zMIC, k):
    term = (a / zMIC)**k
    return psi_max - ((psi_max - psi_min) * term) / (term + 1)

# function in which growht rate is calculated depending on the environmental conditions
def dX_dt(X, t, psi_max, psi_min, zMIC, k):
    if t > 5:
        a_t = 3
    else:
        a_t = 0
    # E is the environmental variation at time t
    current_env = E[int(t) % T]
    growth_rate_modifier = psi_max * reaction_norm(0.3, 0.5, current_env)
    growth_rate = np.log(10) * psi(a_t, growth_rate_modifier, psi_min, zMIC, k) * X

    return max(growth_rate, -X / 0.04)


# Model parameters
psi_min = -2
zMIC = 2
k = 0.8  # Using a single mean k value
psi_max = 0.8  # Assuming this is the psi_max for the simulation

# Time vector
t = np.linspace(0, 10, 10)
# Initial conditions
X0 = 1e3
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
plt.savefig('simple_mod_3.png')
