import numpy as np
# ╔══════════════════════════════════════════════════╗
# ║                  Equations                       ║
# ╚══════════════════════════════════════════════════╝

def environmental_variation(A, B, t, L, R, epsilon):
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
    return (psi_max - psi_min) * (term / (term - psi_min/psi_max))
    return psi_max - ((psi_max - psi_min) * term) / (term + 1) # Giorgio's Implementation

def is_time_for_administration(time):
    # not statement reverses the antibiotic exposure time frames (simply put in front of expression)
    return not time % 10 < 5

def population_is_below_threshold(X, threshold):
    return X < threshold

def growth_rate_modifier(psi_max, params, env):
    return psi_max * reaction_norm(params["I0"], params["b"], env)

def death_rate_modifier(growth):
    return  - growth * 1.5

def realized_variation_calculator(env,X):
    return env
    return env * (1-X/1e9)

def dX_dt(X, t, psi_max, psi_min, zMIC, k, params, environment,antibody_concentration):
    '''function in which growth rate is calculated depending on the environmental conditions'''

    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(t): 
        a_t = antibody_concentration 
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # Environmental variation (as an environmental Cue) at time t
    modified_current_env = realized_variation_calculator(current_env,X)

    modified_growth_rate = growth_rate_modifier(psi_max, params, modified_current_env)
    modified_death_rate = death_rate_modifier(modified_growth_rate)
    growth_rate_after_antibiotic = modified_growth_rate -  psi(a_t, modified_growth_rate, modified_death_rate, zMIC, k)
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return max(actual_growth_rate, -X / 0.04)

def dENV_dt(variables, t, psi_max, psi_min, zMIC, k, params, environment,antibody_concentration):
    '''function in which growth rate is calculated depending on the environmental conditions'''
    X = variables[0]
    realized_variation = variables[1]

    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(t): 
        a_t = antibody_concentration 
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # Environmental variation (as an environmental Cue) at time t
    modified_variation = realized_variation_calculator(current_env, X)
    realized_variation = modified_variation  - realized_variation # substract previous value  in order to reflect the current

    modified_growth_rate = growth_rate_modifier(psi_max, params, modified_variation)
    modified_death_rate = death_rate_modifier(modified_growth_rate)
    growth_rate_after_antibiotic = modified_growth_rate -  psi(a_t, modified_growth_rate, modified_death_rate, zMIC, k)
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return [max(actual_growth_rate, -X / 0.04), realized_variation ]

def dRESPONSE_dt(variables, t, psi_max, psi_min, zMIC, k, params, environment,antibody_concentration):
    '''function in which growth rate is calculated depending on the environmental conditions'''
    X = variables[0]
    actual_response = variables[1]

    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(t): 
        a_t = antibody_concentration 
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # Environmental variation (as an environmental Cue) at time t
    modified_variation = realized_variation_calculator(current_env, X)


    actual_response = reaction_norm(params["I0"], params["b"], modified_variation) - actual_response

    modified_growth_rate = growth_rate_modifier(psi_max, params, modified_variation)
    modified_death_rate = death_rate_modifier(modified_growth_rate)
    growth_rate_after_antibiotic = modified_growth_rate -  psi(a_t, modified_growth_rate, modified_death_rate, zMIC, k)
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return [max(actual_growth_rate, -X / 0.04), actual_response ]