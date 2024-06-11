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

def is_time_for_mutation(time, timepoint):
    return timepoint == time 

def population_is_below_threshold(X, threshold):
    return X < threshold

def growth_rate_modifier(psi_max, response):
    return psi_max * response

def death_rate_modifier(psi_max):
    return  - psi_max  * 1.5

def realized_variation_calculator(env,X):
    return env
    return env * (1-X/1e9)

def sim(initial_conditions, time, env_params, gene_params, antibiotic_framework_params):

    X = initial_conditions[0]

    psi_max = antibiotic_framework_params["psi max"]
    psi_min = antibiotic_framework_params["psi min"]
    antibody_concentration = antibiotic_framework_params["Antibiotic Concentration"]
    zMIC = antibiotic_framework_params["zMIC"]
    k = antibiotic_framework_params["k"]
    A, B, L, R = env_params
    variation_max = A
    variation_min = - variation_max


    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(time): 
        a_t = antibody_concentration 
    else:
        a_t = 0

    

    true_env_variation = environmental_variation(env_params, time) # env variation at time t

    normalized_variation = (true_env_variation - variation_min) / (variation_max - variation_min) 

    theoritical_response = reaction_norm(gene_params, normalized_variation) # response based on variation

    

    modified_psi_max = growth_rate_modifier(psi_max, theoritical_response) # effect of response to psiMax
    modified_death_rate = death_rate_modifier(modified_psi_max) # effect of new psiMax to psiMin
    antibiotic_effect = psi(a_t, modified_psi_max, modified_death_rate, zMIC, k) # effect of antibiotic based on new psiMax/Min

    growth_rate_after_antibiotic = modified_psi_max - antibiotic_effect
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return [max(actual_growth_rate, -X / 0.04), 
            true_env_variation, 
            theoritical_response,
            modified_psi_max, 
            modified_death_rate,
            growth_rate_after_antibiotic,
            antibiotic_effect,
            normalized_variation
            ]

def sim_mutation(initial_conditions, time, env_params, gene_params, antibiotic_framework_params):

    wild_pop = initial_conditions[0]
    mutant_pop = initial_conditions[1]


    psi_max = antibiotic_framework_params["psi max"]
    psi_min = antibiotic_framework_params["psi min"]
    antibody_concentration = antibiotic_framework_params["Antibiotic Concentration"]
    zMIC = antibiotic_framework_params["zMIC"]
    k = antibiotic_framework_params["k"]
    A, B, L, R = env_params
    variation_max = A
    variation_min = - variation_max
    mutated_zMIC = zMIC * 10

    
    if population_is_below_threshold(wild_pop,10):
        wild_pop = 0

    # if population_is_below_threshold(mutant_pop,10):
    #     mutant_pop = 0



    if is_time_for_mutation(int(time),20):
        # print("mutation happened")
        # print(mutation_happened)
        mutant_pop = 10


    if is_time_for_administration(time): 
        a_t = antibody_concentration 
    else:
        a_t = 0    

    true_env_variation = environmental_variation(env_params, time) # env variation at time t
    normalized_variation = (true_env_variation - variation_min) / (variation_max - variation_min) 
    theoritical_response = reaction_norm(gene_params, normalized_variation) # response based on variation    

    modified_psi_max = growth_rate_modifier(psi_max, theoritical_response) # effect of response to psiMax
    modified_death_rate = death_rate_modifier(modified_psi_max) # effect of new psiMax to psiMin

    # ANTIBIOTIC EFFECT ON WILD & MUTANT POPULATIONS
    wild_antibiotic_effect = psi(a_t, modified_psi_max, modified_death_rate, zMIC, k) # effect of antibiotic based on new psiMax/Min
    mutant_antibiotic_effect = psi(a_t, modified_psi_max, modified_death_rate, mutated_zMIC, k) # effect of antibiotic based on new psiMax/Min

    # GROWTH RATES FOR WILD & MUTANT POPULATIONS
    wild_growth_rate_after_antibiotic = modified_psi_max - wild_antibiotic_effect
    mutant_growth_rate_after_antibiotic = modified_psi_max - mutant_antibiotic_effect
    
    # ACTUAL GROWTH FOR WILD & MUTANT POPULATIONS
    wild_actual_growth_rate = np.log(10) * wild_growth_rate_after_antibiotic * wild_pop * (1 - (wild_pop/1e9))
    mutant_actual_growth_rate = np.log(10) * mutant_growth_rate_after_antibiotic * mutant_pop * (1 - (mutant_pop/1e9))


    return [max(wild_actual_growth_rate, -wild_pop / 0.04),
            max(mutant_actual_growth_rate, -mutant_pop / 0.04),
            true_env_variation, 
            theoritical_response,
            modified_psi_max, 
            modified_death_rate,
            wild_growth_rate_after_antibiotic,
            wild_antibiotic_effect,
            normalized_variation
            ]
