# ╔══════════════════════════════════════════════════╗
# ║                  Parameters                      ║
# ╚══════════════════════════════════════════════════╝
import numpy as np
from tools import *

#region

# All environments must have different keys otherwise will be overwritten
# All environments must have at least one different value otherwise only the last will be saved
determistic = [0.3,0.6,0.9]
stochastic = [0.0,]
lifespan = [10]
relativeVariation = [1,3]
timesteps = [101]

environments_params = construct_params(determistic, stochastic, lifespan, relativeVariation, timesteps)


genotypes_params = {
    "Genotype 1": {"I0": 0.1, "b": 0.9},
    "Genotype 3": {"I0": 0.5, "b": 0.4},
    "Genotype 4": {"I0": 0.8, "b": 0},    
}

# antibiotic response curve
antibody_concentration = 100
psi_min = -2 # maximum death rate
zMIC = 2 # concentration in which net growth rate is zero
k = 0.8  # Using a single mean k value
psi_max = 0.3  # maximal growth rate
initial_populations = [1e7]



# for all simulations and layer appliers to work properly
# the slicing must be at least time+1 (e.g. 101 slices for time=100)
time_frame = np.linspace(0, 100, 101) #should be passed on odeint()


#endregion