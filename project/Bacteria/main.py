import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.gridspec import GridSpec

def psi(a, psi_max, psi_min, zMIC, kappa, k):
    term = (a / zMIC)**k
    return psi_max - ((psi_max - psi_min) * term) / (term - (psi_min / psi_max))

def a(t):
    # antibiotic decay function with half-life of 1 hour
    if t > 6 :
        return 12.5
    else:   
        return 0 

def dX_dt(X, t, psi_max, psi_min, zMIC, kappa, k):
    a_t = a(t)
    growth_rate = np.log(10) * psi(a_t, psi_max, psi_min, zMIC, kappa, k) * X
    return max(growth_rate, -X / 0.04)

## renders list of different dose concentration responses a(c)
def dose_response(a_values, psi_max, psi_min, zMIC, kappa, k):
    responses = np.zeros(len(a_values))
    for i, a_val in enumerate(a_values):
        responses[i] = psi(a_val, psi_max, psi_min, zMIC, kappa, k)
    return responses

# runs dx/dt for numRuns
def run_simulation(mean_k, std_k, psi_max_mutant, ax_dynamics, ax_dose, ax_summed, label, color, num_runs=100):
    final_densities = np.zeros(num_runs) # renders zeros []
    k_values = np.random.lognormal(mean=mean_k, sigma=std_k, size=num_runs) # [different steepnesses of antibiotic response curve]
    a_values = np.linspace(0, 40, 100) #  [different concentrations of antibiotic]
    all_densities = np.zeros((num_runs, len(t))) # (arg is tuple (int,int)) renders an [runs x timesteps] array

    # different starting population for mutants
    if psi_max_mutant == 0.7:
        X0_mutant = 1e3
    else:
        X0_mutant = 1e6

    # grab k and solve dx/dt for this k for numRuns
    for run in range(num_runs):
        k = k_values[run]
        X = odeint(dX_dt, X0_mutant, t, args=(psi_max_mutant, psi_min, zMIC_mutant, kappa, k)) # dx/dt
        final_density = max(X[-1, 0], 0) #X[] expands the list here??
        final_densities[run] = final_density
        all_densities[run, :] = X[:, 0]

        if run < 20:  # Plot only the first 20 runs to avoid clutter
            ax_dynamics.plot(t, X, color=color, alpha=0.7, linewidth=1)
            responses = dose_response(a_values, psi_max_mutant, psi_min, zMIC_mutant, kappa, k)
            ax_dose.plot(a_values, responses, color=color, alpha=0.7, linewidth=1)
            
    # Plot summed densities
    summed_densities = np.sum(all_densities, axis=0) / num_runs # mean densities here??
    ax_summed.plot(t, summed_densities, label=label, color=color, linewidth=2)

    # why call them again??
    ax_dynamics.plot([], [], color="black", label=label, linewidth=2)
    ax_dose.plot([], [], color=color, label=label, linewidth=2)
    return final_densities, k_values, summed_densities

# Model parameters
psi_min = -2
kappa = 1.0
t = np.linspace(0, 10, 100)

# Simulation parameters
mean_k_values = [0.8]
std_values = [0.1, 0.8]
zMIC_values = [4, 8]
colors = ['blue', 'green']
psi_max_values = [0.8, 0.7]

for zMIC_mutant in zMIC_values:
    for mean_k in mean_k_values:
        # Create plots
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(2, 2, height_ratios=[1, 1])
        ax_dynamics = fig.add_subplot(gs[0, 0])
        ax_dynamics.set_xlabel('Time')
        ax_dynamics.set_ylabel('Bacterial Density (Log Scale)')
        ax_dynamics.set_yscale('log')
        ax_dynamics.set_ylim(1, 1e9)
        ax_dynamics.set_title(f'Bacterial Growth Dynamics (zMIC={zMIC_mutant}, mean_k={mean_k})')

        ax_dose = fig.add_subplot(gs[0, 1])
        ax_dose.set_xlabel('Antibiotic Concentration (a)')
        ax_dose.set_ylabel('Growth Rate Factor (ψ)')
        ax_dose.set_title(f'Dose-Response Curves (zMIC={zMIC_mutant}, mean_k={mean_k})')
        ax_dose.axvline(x=zMIC_mutant, color='magenta', linestyle='--', label='MIC')

        ax_summed = fig.add_subplot(gs[1, 0])
        ax_summed.set_xlabel('Time')
        ax_summed.set_ylabel('Summed Bacterial Density')
        ax_summed.set_yscale('log')
        ax_summed.set_ylim(1, 1e9)
        ax_summed.set_title(f'Summed Bacterial Densities Over Time (zMIC={zMIC_mutant}, mean_k={mean_k})')

        ax_k_distribution = fig.add_subplot(gs[1, 1])
        ax_k_distribution.set_title(f'Distribution of k (zMIC={zMIC_mutant}, mean_k={mean_k})')
        ax_k_distribution.set_xlabel('k Value')
        ax_k_distribution.set_ylabel('Probability Density')
        ax_k_distribution.set_xlim(0, 40)

        all_k_values = []
        all_final_densities = []
        all_summed_densities = []

        for i, (std_k, psi_max_mutant, color) in enumerate(zip(std_values, psi_max_values, colors)):
            label = f'Std={std_k:.1f}, Psi_max={psi_max_mutant}'
            final_densities, k_values, summed_densities = run_simulation(mean_k, std_k, psi_max_mutant, ax_dynamics, ax_dose, ax_summed, label, color)
            all_k_values.append(k_values)
            all_final_densities.append(final_densities)
            all_summed_densities.append(summed_densities)

        # Plot the sum of bacterial densities for both subpopulations as a dotted line
        total_summed_densities = np.sum(all_summed_densities, axis=0)
        ax_summed.plot(t, total_summed_densities, 'k--', label='Total Summed Densities')

        # Plot histograms
        for i, (std_k, k_values, final_densities, color) in enumerate(zip(std_values, all_k_values, all_final_densities, colors)):
            final_densities[final_densities < 1] = 0
            num_bins = 100
            bins = np.logspace(0, 10, num_bins)
            ax_k_distribution.hist(k_values, bins=num_bins, density=True, alpha=0.5, color=color, label=f'Std={std_k:.1f}')
            #ax_histogram.hist(final_densities, bins=bins, density=True, alpha=0.5, color=color, label=f'Std={std_k:.1f}')

        ax_dynamics.legend()
        ax_dose.legend()
        ax_summed.legend()
        ax_k_distribution.legend()
        #ax_histogram.legend()
        plt.tight_layout()
        plt.savefig(f'Example_zMIC{zMIC_mutant}_meanK{mean_k}.png')