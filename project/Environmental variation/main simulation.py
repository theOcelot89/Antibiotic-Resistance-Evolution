from classes import *
# ╔══════════════════════════════════════════════════╗
# ║                  Simulations                     ║
# ╚══════════════════════════════════════════════════╝

# region test simulations
environment = Environment()
environment.trim()
environment.save()

environment.gene_reaction_norms(genotypes_params)
environment.gene_responses(genotypes_params)
environment.run_simulation(genotypes_params, initial_populations)

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







