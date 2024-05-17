from utils.classes import *
from utils.parameters import *


# ╔══════════════════════════════════════════════════╗
# ║                  Simulations                     ║
# ╚══════════════════════════════════════════════════╝
A, B, L, R, t = 1.0, 0, 10, 8, 100

# region environment simulations
environment = Environment(A , B , L , R , t , genotypes_params, antibiotic_framework)
environment.variation()
environment.gene_reaction_norms()
environment.responses()
environment.dynamics()
environment.dynamics_with_antibiotic_frames()
environment.dynamics_with_antibiotic_frames_and_variation()
environment.actual_psi_max()
environment.growth_rate_after_antibiotic()
environment.actual_psi_max__antibiotic_effect__growth_rate()
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







