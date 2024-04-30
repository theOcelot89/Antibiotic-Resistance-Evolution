from utils.classes import *
from utils.parameters import *


# ╔══════════════════════════════════════════════════╗
# ║                  Simulations                     ║
# ╚══════════════════════════════════════════════════╝

# region environment simulations
environment = Environment(A = 0.9, B = 0, L = 10, R = 8, t = 200)
environment.trim()
environment.save()
# environment.gene_reaction_norms(genotypes_params)
# environment.gene_responses(genotypes_params)
# environment.population_dynamics(genotypes_params,antibiotic_framework)
# environment.population_dynamics_antibiotic_frames(genotypes_params,antibiotic_framework)
environment.population_dynamics_antibiotic_frames_env_variation(genotypes_params, antibiotic_framework)
# environment.run_simulation(genotypes_params, antibiotic_framework)
environment.realized_variation(genotypes_params, antibiotic_framework)
environment.actual_response(genotypes_params,antibiotic_framework)
#endregion

#region main simulations
# simulator = Simulator(environments_params, genotypes_params, antibiotic_framework)
# simulator.yield_environment_plots()
# simulator.yield_phenotypic_responses()
# simulator.yield_reaction_norms()
# simulator.yield_population_dynamics()
# simulator.yield_environment_plots_with_antibiotic_frames()
# simulator.yield_population_dynamics_with_antibiotic_frames()
# simulator.yield_population_dynamics_with_antibiotic_frames_env_variation()
# simulator.run()
#endregion







