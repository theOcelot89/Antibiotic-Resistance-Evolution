from parameters import *
from tools import *
from equations import *
from classes import *




# ╔══════════════════════════════════════════════════╗
# ║                  Simulations                     ║
# ╚══════════════════════════════════════════════════╝

# region test simulations

    #region environment construction
environment = Environment()
environment.trim()
environment.save()
    #endregion

    #region norms & responses to environmental variation
environment.gene_reaction_norms(genotypes_params)
environment.gene_responses(genotypes_params)
    #endregion

    #region bacterial growth simulations
environment.run_simulation(genotypes_params, initial_populations)
    #endregion

#endregion

#region main simulations
simulator = Simulator(environments_params, genotypes_params)
simulator.yield_environment_plots()
simulator.yield_phenotypic_responses()
simulator.yield_reaction_norms()
simulator.yield_population_dynamics()
simulator.yield_environment_plots_with_antibiotic_frames()
simulator.yield_population_dynamics_with_antibiotic_frames()
simulator.yield_population_dynamics_with_antibiotic_frames_env_variation()
# simulator.generate_report()
# simulator.run()
#endregion







