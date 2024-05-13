import numpy as np
from .tools import construct_params
# ╔══════════════════════════════════════════════════╗
# ║                  Parameters                      ║
# ╚══════════════════════════════════════════════════╝
# All environments must have different keys otherwise will be overwritten
# All environments must have at least one different value otherwise only the last will be saved
determistic = [0.4,0.7,1]
stochastic = [0.0]
lifespan = [10]
relativeVariation = [1,2,4]
timesteps = [101]


genotypes_params = {
    "Genotype 1": {"I0": 0.1, "b": 0.9},
    "Genotype 3": {"I0": 0.5, "b": 0.4},
    "Genotype 4": {"I0": 0.8, "b": 0},    
}


#################
#   IMPORTANT!  #
#################
# for all simulations and layer appliers to work properly
# the slicing must be at least time+1 (e.g. np.linspace(0, 100, 101))


antibiotic_framework = {
    "zMIC" : 4,
    "Antibiotic Concentration" : 30,
    "psi_max" : 0.3,
    "psi_min" : -2,
    "k" : 0.8,
    "time frame" : np.linspace(0, 100, 101),
    "Initial Populations": [1e7]
}



environments_params = construct_params(determistic, stochastic, lifespan, relativeVariation, timesteps)


if __name__ == "__main__":

    for name, parameters in environments_params.items():
        print(name, parameters)

    for name, parameters in genotypes_params.items():    
        print(name, parameters)