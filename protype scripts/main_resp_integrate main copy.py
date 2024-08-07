# ╔══════════════════════════════════════════════════╗
# ║                   Imports                        ║
# ╚══════════════════════════════════════════════════╝
#region
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                   Classes                        ║
# ╚══════════════════════════════════════════════════╝
#region

class Environment():

    def __init__(self, A = 1, B = 0.1, L = 10, R = 10, t = 365):
        '''
        Environmental variation that individuals face, relative to their lifespam
        t = time, A = determinism magnitude, B = stochasticity magnitude
        L = lifespan, R = generations/environmental cycle, epsilon = stochastic error term
        https://doi.org/10.1073/pnas.1408589111
        '''
        self.A = A 
        self.B = B
        self.t = np.arange(t)
        self.L = L
        self.R = R
        self.epsilon = np.random.normal(0, 1, t)
        self.trimmed = False # flag for trimming to put in the plot's title

        # yield variation
        self.variation = self.A * np.sin(2 * np.pi * self.t / (self.L * self.R)) + self.B * self.epsilon
        # construct plot and immediately unpack
        self.fig, self.ax = self._create_plot()    
        
        
        print(f"New environment created!")

    def _create_plot(self):
            
        fig, ax = plt.subplots(figsize=(12, 6))   
        ax.plot(self.t, self.variation, label='Environmental Variation')
        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 
        ax.set_title(f' Environmental variation A={self.A}, B={self.B}, L={self.L}, R={self.R}, trimmed={self.trimmed} ')

        return fig, ax

    def trim(self):
        '''limiting environmental variation between 0 & 1'''
        # because trim will work directly on self.variation its action is irreversible & will follow the instance for the rest of the script when plots are constructed 

        self.variation[self.variation < 0] = 0  # if any of the values are negative, set them to 0
        self.variation[self.variation > 1] = 1    # if any of the values are greater than 1, set them to 1
        self.trimmed = True

        # reconstruct variation plot
        self.fig, self.ax = self._create_plot()

    def view(self):

        renderPeriod = 3
        initial_title = self.ax.get_title() # keep original title to put after rendering
        self.ax.set_title(self.ax.get_title() + "\n" + f'Will disappear after {renderPeriod} seconds') # temp title for rendering
        plt.show(block=False) # stop from blocking in the execution
        plt.pause(renderPeriod)
        self.ax.set_title(initial_title) # set initial title again
        
    def save(self):
        self.fig.savefig(f"Env. variation t={len(self.t)}, A={self.A}, B={self.B}, L={self.L}, R={self.R}, trimmed={self.trimmed}.png")

    def gene_responses(self, genotypes):

        fig, ax = self._create_plot() # create a copy of the current variation plot (keep clean the original)

        for name, params in genotypes.items():
            I = reaction_norm(params["I0"], params["b"], self.variation)
            ax.plot(self.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")

        ax.set_title('Phenotypic Response')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Phenotypic response (I)')
        ax.legend(bbox_to_anchor=(1.34, 1))

        fig.savefig("Responses to Environmental variation")

    def gene_reaction_norms(self, genotypes):

        fig, ax = plt.subplots(figsize=(12,6)) # create plot from scratch

        for name, params in genotypes.items():
            I = reaction_norm(params["I0"], params["b"], self.variation)
            ax.plot(self.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")

        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot 

        ax.set_title('Reaction Norms')
        ax.set_xlabel('Environmental Variation (E)')
        ax.set_ylabel('Phenotypic response (I)')

        fig.savefig("Reaction Norms")

    def run_simulation(self, genotypes, initial_populations):

        fig, ax = plt.subplots(figsize=(14,6)) # prepare plot

        for initial_population in initial_populations:
            for name, params in genotypes.items():

                X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, self)) # args will be passed down to dX_dt
                ax.plot(time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')

        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e9)   

        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.41, 1), fontsize="7") # place legend out of plot

        fig.savefig(f'Genotypes dynamics.png')

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Functions                       ║
# ╚══════════════════════════════════════════════════╝
#region

def yield_environments(environments_parameters): 

    environment_list = []

    for name, params in environments_parameters.items():
        A, B, L, R, t = params.values() # unpacking env parameters
        env = Environment(A, B, L, R, t)
        env.trim()
        env.save()
        environment_list.append(env)
    
    return environment_list

def yield_environment_plots(environments):
    # this technique is based on matplots basic tutorial
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html

    fig = plt.figure(figsize=(12, len(environments)*5)) # empty figure for template
    gs = fig.add_gridspec(len(environments), hspace=0) # grid with dimensions & space between plots
    axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
    fig.suptitle(u"Environmental Variations E\u209C = A·sin(2πt/LR) + B·ε", fontsize = 30)
    for ax in axs.flat:
        ax.set(xlabel='Time (t)', ylabel='Environmental variation (E)')

    for i, env in enumerate(environments):
        axs[i].plot(env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}')
        axs[i].legend()
        # axs[i].set_ylim(0,1)

    fig.savefig("Stacked Environments")

def yield_phenotypic_responses(environments, genotypes):

    fig = plt.figure(figsize=(12, len(environments)*5)) # empty figure for template, dynamic height of plot
    gs = fig.add_gridspec(len(environments), hspace=0) # grid with dimensions & space between plots
    axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
    fig.suptitle(u"Phenotypic Responses To  Environmental Variations", fontsize = 30)
    for ax in axs.flat:
        ax.set(xlabel='Time (t)', ylabel='Response (I)')

    for i, env in enumerate(environments):

        for name, params in genotypes.items():
            I = reaction_norm(params["I0"], params["b"], env.variation)
            axs[i].plot(env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
            axs[i].legend()

    fig.savefig("Stacked Phenotypic Responses")

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
    return psi_max - ((psi_max - psi_min) * term) / (term + 1)

def dX_dt(X, t, psi_max, psi_min, zMIC, k, params, environment):
    '''function in which growth rate is calculated depending on the environmental conditions'''
    t_g = 18
    ant_a = 5
    # decide in which timestep(e.g day) to quit the administration of antibiotic
    if t % 18 < 5:
        a_t = 2.3 # antibiotic concentration
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # Environmental variation (as an environmental Cue) at time t
    growth_rate_modifier = psi_max * reaction_norm(params["I0"], params["b"], current_env) # new psimax depending on plasticity
    deathrateModifier = - (growth_rate_modifier * 5)
    growth_rate = np.log(10) * psi(a_t, growth_rate_modifier, deathrateModifier, zMIC, k) * X * (1 - (X/1e9))

    return max(growth_rate, -X / 0.04)

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Parameters                      ║
# ╚══════════════════════════════════════════════════╝
#region Environment
# All environments must have different keys otherwise will be overwritten
# All environments must have at least one different value otherwise only the last will be saved
environments_params = {
    "Env 1": {"A": 0.3, "B": 0.0, "L": 10, "R": 2, "t": 110},
    "Env 2": {"A": 0.6, "B": 0.0, "L": 10, "R": 2, "t": 110},
    "Env 3": {"A": 0.9, "B": 0.0, "L": 10, "R": 2, "t": 110},
}
#endregion

#region Genotypes
genotypes_params = {
    "Genotype 1": {"I0": 0.25, "b": 1},
    "Genotype 2": {"I0": 0.3, "b":0.5},
    "Genotype 3": {"I0": 0.6, "b": 0.05},
}
#endregion

#region Antibiotic Response Curve

psi_min = -2 # maximum death rate
zMIC = 2 # concentration in which net growth rate is zero
k = 0.8  # Using a single mean k value
psi_max = 0.8  # maximal growth rate

time_frame = np.linspace(0, 50, 50) #should be passed on odeint()
initial_populations = [1e7]

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Simulations                     ║
# ╚══════════════════════════════════════════════════╝

#region test simulations


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
environments = yield_environments(environments_params)  # create several environments
yield_environment_plots(environments) # multiplot for all environments
yield_phenotypic_responses(environments, genotypes_params) # multiplot for phenotypic respones

#endregion

#region population growth multiplot
fig = plt.figure(figsize=(12, len(environments)*5)) # empty figure for template, dynamic height of plot
gs = fig.add_gridspec(len(environments), hspace=0) # grid with dimensions & space between plots
axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
fig.suptitle(f"Population Dynamics \nResponse Curve Parameters: k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}", fontsize = 20)
for ax in axs.flat:
    ax.set(xlabel='Time (t)', ylabel='Bacterial Density')

for i, env in enumerate(environments):

    for initial_population in initial_populations:
        for name, params in genotypes_params.items():

            X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env)) # args will be passed down to dX_dt
            axs[i].plot(time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}')
            axs[i].set_yscale('log')
            axs[i].set_ylim(1, 1e9) 
            axs[i].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")  
fig.savefig("Stacked Population Dynamics")

#endregion