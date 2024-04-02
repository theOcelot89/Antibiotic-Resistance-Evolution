# ╔══════════════════════════════════════════════════╗
# ║                  Imports                         ║
# ╚══════════════════════════════════════════════════╝
#region
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from typing import Annotated, Any
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Functions                       ║
# ╚══════════════════════════════════════════════════╝
#region


def environmental_variation(t, A, B, L, R, epsilon):
    '''
    Environmental variation that individuals face, relative to their lifespam
    t = time, A = determinism magnitude, B = stochasticity magnitude
    L = lifespan, R = generations/environmental cycle, epsilon = stochastic error term
    https://doi.org/10.1073/pnas.1408589111
    '''

    return A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon

def reaction_norm(I0, b, C):
    '''
    Estimation of individuals' phenotypic reaction to environmental variation 
    I0 = baseline amount, b = degree of plasticity, C = environmental cue
    https://doi.org/10.1073/pnas.1408589111
    '''

    return I0 + b * C

def psi(a, psi_max, psi_min, zMIC, k):
    '''
    Effect of antibiotic on bacterial growth
    a = concentration of antibiotic, psiMax = max growth rate, psiMin = max decline rate
    zMIC = Minimun Inhibitory Concentration, k = steepness of the antibiotic response curve
    https://doi.org/10.1371/journal.pcbi.1011364
    '''

    term = (a / zMIC)**k
    return psi_max - ((psi_max - psi_min) * term) / (term + 1)

def dX_dt(X, t, psi_max, psi_min, zMIC, k):
    '''function in which growth rate is calculated depending on the environmental conditions'''

    # decide in which timestep(e.g day) to quit the administration of antibiotic
    if t > 5: 
        a_t = 3 # antibiotic concentration
    else:
        a_t = 0
    
    current_env = E[int(t) % T] # E is the environmental variation (as an environmental Cue) at time t
    growth_rate_modifier = psi_max * reaction_norm(params["I0"], params["b"], current_env) # new psimax depending on plasticity
    growth_rate = np.log(10) * psi(a_t, growth_rate_modifier, psi_min, zMIC, k) * X

    return max(growth_rate, -X / 0.04)

def E_fixation(env):
    '''limiting environmental variation between 0 & 1'''

    env[env < 0] = 0  # if any of the values are negative, set them to 0
    env[env > 1] = 1    # if any of the values are greater than 1, set them to 1

    return E


#endregion

# ╔══════════════════════════════════════════════════╗
# ║                 Parameters                       ║
# ╚══════════════════════════════════════════════════╝
#region Environmental variation


T = 110 # number of time steps to simulate (e.g., days in a year)
epsilon = np.random.normal(0, 1, T) # stochastic error term
A = 1 # determinism magnitude
B = [0.1, 0.3, 0.9] # stochasticity magnitude
L = 10 # lifespan
R =  100 # generations/environmental cycle




    

#endregion

#region Genotypes

genotypes = {
    "Genotype 1": {"I0": 0, "b": 0},
    "Genotype 2": {"I0": 0, "b":0.5},
    "Genotype 3": {"I0": 0, "b": 1}
}
colors = ['blue', 'green', 'red'] # Define colors for each genotype

#endregion

#region antibiotic response curve

psi_min = -2 # maximum death rate
zMIC = 2 # concentration in which net growth rate is zero
k = 0.8  # Using a single mean k value
psi_max = 0.8  # maximal growth rate

# Time vector
t = np.linspace(0, 10, 10)
initial_populations = [1e3,1e4,1e5,1e6]

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                 Plots/Simulations                ║
# ╚══════════════════════════════════════════════════╝

#region responses to environmental variation

for B in B: # different stochasticity amplifiers

    E = environmental_variation(np.arange(T), A, B, L, R, epsilon)

    E = E_fixation(E)

    fig, ax = plt.subplots(figsize=(12, 6))   

    # Plot environmental variation in time
    ax.plot(np.arange(T), E, label='Environmental Variation', linestyle='--', color='grey')
    pos = ax.get_position() #returns bbox in order to manipulate width/height
    ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
    ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 

    fig.savefig(f'b = {B}environmental_variation.png')

    # plot individuals' response to environmental variation
    for (name, params), color in zip(genotypes.items(), colors):
        I = reaction_norm(params["I0"], params["b"], E)
        ax.plot(np.arange(T), I, label=f"{name},IO={params["I0"]},b ={params["b"]}", color=color)

    ax.set_title('Phenotypic Response')
    ax.set_xlabel('Time (t)')
    ax.set_ylabel('Phenotypic response (I)')
    ax.legend(bbox_to_anchor=(1.32, 1))

    # Save the plot to a file
    fig.savefig(f'for{B} stochasticity Reaction Norms Across Time.png')

#endregion

#region reaction norms
    
fig, ax = plt.subplots(figsize=(12,6))

for (name, params), color in zip(genotypes.items(), colors):
    I = reaction_norm(params["I0"], params["b"], E)
    plt.plot(E, I, label=f"{name},IO={params["I0"]},b ={params["b"]}", color=color)

pos = ax.get_position() #returns bbox in order to manipulate width/height
ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 

ax.set_title('Reaction Norms for Genotypes')
ax.set_xlabel('Environmental Variation (E)')
ax.set_ylabel('Phenotypic Expression (I)')
# ax.legend()

fig.savefig('Reaction Norms for Genotypes.png') # Save the plot


#endregion

#region bacterial growth simulations

    #region different initial populations

fig, ax = plt.subplots(figsize=(14,6))

for X0 in initial_populations:

    X = odeint(dX_dt, X0, t, args=(psi_max, psi_min, zMIC, k)) # args will be passed down to dX_dt
    
    ax.plot(t, X, label=f'X0={'{:.0e}'.format(X0)} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')
    ax.set_xlabel('Time')
    ax.set_ylabel('Bacterial Density')
    ax.set_yscale('log')
    ax.set_ylim(1, 1e9)

pos = ax.get_position() #returns bbox in order to manipulate width/height
ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
ax.legend(bbox_to_anchor=(1.4, 1), fontsize="7") # place legend out of plot

fig.savefig(f' Different initial population dynamics.png')

    #endregion

    #region different genotypes

fig, ax = plt.subplots(figsize=(14,6))

for (name, params), color in zip(genotypes.items(), colors):
    
    X = odeint(dX_dt, initial_populations[0], t, args=(psi_max, psi_min, zMIC, k)) # args will be passed down to dX_dt
    
    ax.plot(t, X, label=f'X0={'{:.0e}'.format(initial_populations[0])} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')
    ax.set_xlabel('Time')
    ax.set_ylabel('Bacterial Density')
    ax.set_yscale('log')
    ax.set_ylim(1, 1e9)

pos = ax.get_position() #returns bbox in order to manipulate width/height
ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
ax.legend(bbox_to_anchor=(1.4, 1), fontsize="7") # place legend out of plot

fig.savefig(f' Different genotypes dynamics.png')
    #endregion

#endregion