import numpy as np
# ╔══════════════════════════════════════════════════╗
# ║                  Equations                       ║
# ╚══════════════════════════════════════════════════╝

def environmental_variation(params,t):

    A, B, L, R = params
    epsilon = np.random.normal(0, 1)

    return A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon


def reaction_norm(params, C):
    '''
    Estimation of individuals' phenotypic reaction to environmental variation 
    I0 = baseline amount, b = degree of plasticity, C = environmental cue
    https://doi.org/10.1073/pnas.1408589111
    '''

    I0 = params["I0"]
    b = params["b"]

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
    # not statement reverses the antibiotic exposure time frames (simply put "not" in front of expression)
    return not time % 10 < 5

def population_is_below_threshold(X, threshold):
    return X < threshold

def growth_rate_modifier(psi_max, response):
    return psi_max * response

def death_rate_modifier(psi_max):
    return  - psi_max  * 1.5

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
    
    A = environment.A
    B = environment.B
    L = environment.L
    R = environment.R


    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(t): 
        a_t = antibody_concentration 
    else:
        a_t = 0
    
    epsilon = np.random.normal(0, 1)
    true_env_variation = environmental_variation(A, B, t, L, R, epsilon)

    theoritical_response = reaction_norm(params["I0"], params["b"], true_env_variation)

    modified_growth_rate = growth_rate_modifier(psi_max, params, true_env_variation)
    modified_death_rate = death_rate_modifier(modified_growth_rate)
    growth_rate_after_antibiotic = modified_growth_rate -  psi(a_t, modified_growth_rate, modified_death_rate, zMIC, k)
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return [max(actual_growth_rate, -X / 0.04), true_env_variation, theoritical_response]

def sim(initial_conditions, time, env_params, gene_params, antibiotic_framework_params):

    X = initial_conditions[0]

    psi_max = antibiotic_framework_params["psi max"]
    psi_min = antibiotic_framework_params["psi min"]
    antibody_concentration = antibiotic_framework_params["Antibiotic Concentration"]
    zMIC = antibiotic_framework_params["zMIC"]
    k = antibiotic_framework_params["k"]

    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(time): 
        a_t = antibody_concentration 
    else:
        a_t = 0

    true_env_variation = environmental_variation(env_params, time) # env variation at time t
    theoritical_response = reaction_norm(gene_params, true_env_variation) # response based on variation

    actual_psi_max = growth_rate_modifier(psi_max, theoritical_response) # effect of response to psiMax
    modified_death_rate = death_rate_modifier(actual_psi_max) # effect of new psiMax to psiMin
    antibiotic_effect = psi(a_t, actual_psi_max, modified_death_rate, zMIC, k) # effect of antibiotic based on new psiMax/Min

    growth_rate_after_antibiotic = actual_psi_max - antibiotic_effect
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return [max(actual_growth_rate, -X / 0.04), 
            true_env_variation, 
            theoritical_response,
            actual_psi_max, 
            modified_death_rate
            ]