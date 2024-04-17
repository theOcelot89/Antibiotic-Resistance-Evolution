def psi(a, psi_max, psi_min, zMIC, k):
    '''
    Effect of antibiotic on bacterial growth
    a = concentration of antibiotic, psiMax = max growth rate, psiMin = max decline rate
    zMIC = Minimun Inhibitory Concentration, k = steepness of the antibiotic response curve
    https://doi.org/10.1371/journal.pcbi.1011364
    '''

    term = (a / zMIC)**k
    # return (psi_max - psi_min) * (term / (term - psi_min/psi_max))

    return psi_max - ((psi_max - psi_min) * term) / (term + 1)

result = psi(a=0, psi_max= 8, psi_min= -2 ,zMIC= 2, k=0.8)

print(result)