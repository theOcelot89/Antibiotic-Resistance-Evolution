# ╔══════════════════════════════════════════════════╗
# ║                   Imports                        ║
# ╚══════════════════════════════════════════════════╝
#region
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import io 
import os
from PIL import Image 
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
        self.variation = environmental_variation(self.A, self.B, self.t, self.L, self.R, self.epsilon)
        # self.variation = self.A * np.sin(2 * np.pi * self.t / (self.L * self.R)) + self.B * self.epsilon
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

class Simulator():
    '''
    Simulator class is designed to process environments & genotypes.
    It produces plots for environmental variation, genotype responses & norms and population dynamics
    based on some predefined equations.
    '''
    def __init__(self, environment_params, genotype_params):

        self.environment_params = environment_params
        self.genotypes = genotype_params

        self.environments = self._yield_environments() # immediatly create environments
    
    def run(self):

        self.envs_plot = self.yield_environment_plots()
        self.norms_plot = self.yield_reaction_norms()
        self.responses_plot = self.yield_phenotypic_responses()
        self.dynamics_plot = self.yield_population_dynamics()

        self._generate_report()
       
    def _yield_environments(self): 

        environment_list = []

        for name, params in self.environment_params.items():
            A, B, L, R, t = params.values() # unpacking env parameters
            env = Environment(A, B, L, R, t) # create Environment Instance
            env.trim()
            # env.save()
            environment_list.append(env)
        
        return environment_list

    def _generate_report(self):
        # ok that was a touch one & a number of steps & guides have to be considered in order for a nice report sheet
        # 1. convert figures to images (that is the only way to combine them) https://www.geeksforgeeks.org/saving-a-plot-as-an-image-in-python/
        # 2. use fig.figAddSubplot technique to create the grid https://www.geeksforgeeks.org/how-to-display-multiple-images-in-one-figure-correctly-in-matplotlib/
        # 3. use plt.imshow to render the image in the figure https://www.geeksforgeeks.org/how-to-display-multiple-images-in-one-figure-correctly-in-matplotlib/
         
        # print(self.norms_plot.get_size_inches())  #dynamically set width/height from the dimensions of the plots

        img1 = fig2img(self.dynamics_plot)
        img2 = fig2img(self.envs_plot)
        img3 = fig2img(self.responses_plot)
        img4 = fig2img(self.norms_plot)

        fig = plt.figure(figsize=(20,7), layout="compressed")

        fig.add_subplot(141)
        plt.imshow(img1)
        plt.axis('off') 
        fig.add_subplot(142)
        plt.imshow(img2)
        plt.axis('off') 
        fig.add_subplot(143)
        plt.imshow(img3)
        plt.axis('off') 
        fig.add_subplot(144)
        plt.imshow(img4)
        plt.axis('off')

        # plt.subplots_adjust()
        # fig.tight_layout()
        # fig.savefig('Report', dpi=600, bbox_inches='tight') # dpi for a better resolution, bboxinches for trimming margins
        save('./report/report', dpi=600)

    def yield_environment_plots(self):
        # this technique is based on matplots basic tutorial
        # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html

        fig = plt.figure(figsize=(12, len(self.environments)*5)) # empty figure for template
        gs = fig.add_gridspec(len(self.environments), hspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(u"Environmental Variations\nE\u209C = A·sin(2πt/LR) + B·ε", fontsize = 30)
        for ax in axs.flat:
            ax.set(xlabel='Time (t)', ylabel='Environmental variation (E)')

        for i, env in enumerate(self.environments):
            axs[i].plot(env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}')
            axs[i].legend()
            # axs[i].set_ylim(0,1)

        # save('./report/Stacked Environments')
        return fig

    def yield_reaction_norms(self):
        
        fig = plt.figure(figsize=(12, len(self.environments)*5)) # empty figure for template, dynamic height of plot
        gs = fig.add_gridspec(len(self.environments), hspace=0.1) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True,) # sharing the same y range (i think based on the bigger value)
        axs[-1].set_xlabel('Environmental variation (E)') # put label only on the bottom plot
        fig.suptitle(u"Reaction Norms\n I=I\u2080 + b·C", fontsize = 30)
        
        for ax in axs.flat:
            ax.set(ylabel='Response (I)')

        for i, env in enumerate(self.environments):

            for name, params in self.genotypes.items():
                I = reaction_norm(params["I0"], params["b"], env.variation)
                axs[i].plot(env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                axs[i].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
                axs[i].grid(True)

        # save('./report/Stacked Reaction Norms')
        return fig
    
    def yield_phenotypic_responses(self):

        fig = plt.figure(figsize=(12, len(self.environments)*5)) # empty figure for template, dynamic height of plot
        gs = fig.add_gridspec(len(self.environments), hspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(u"Phenotypic Responses\n I\u209C=I\u2080 + b·E\u209C", fontsize = 30)
        for ax in axs.flat:
            ax.set(xlabel='Time (t)', ylabel='Response (I)')

        for i, env in enumerate(self.environments):

            for name, params in self.genotypes.items():
                I = reaction_norm(params["I0"], params["b"], env.variation)
                axs[i].plot(env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                axs[i].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")

        # save('./report/Stacked Phenotypics Responses')
        return fig
    
    def yield_population_dynamics(self):

        fig = plt.figure(figsize=(12, len(self.environments)*5)) # empty figure for template, dynamic height of plot
        gs = fig.add_gridspec(len(self.environments), hspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(f"Population Dynamics \nResponse Curve Parameters: k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}", fontsize = 20)
        for ax in axs.flat:
            ax.set(xlabel='Time (t)', ylabel='Bacterial Density')

        for i, env in enumerate(self.environments):
            for initial_population in initial_populations:
                for name, params in self.genotypes.items():

                    X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env)) # args will be passed down to dX_dt
                    axs[i].plot(time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}')
                    axs[i].set_yscale('log')
                    axs[i].set_ylim(1, 1e9) 
                    axs[i].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")  

        # save('./report/Stacked Population Dynamics')
        return fig
    
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                    Tools                         ║
# ╚══════════════════════════════════════════════════╝
#region

def fig2img(fig): 
    buf = io.BytesIO() 
    fig.savefig(buf) 
    buf.seek(0) 
    img = Image.open(buf) 
    return img 

def save(path, ext='png', close=True, verbose=True, **kwargs):
    # https://gist.github.com/jhamrick/5320734
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.
    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.
    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.
    """
    
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'

    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)

    # The final path to save to
    savepath = os.path.join(directory, filename)

    if verbose:
        print("Saving figure to '%s'..." % savepath),

    # Actually save the figure
    plt.savefig(savepath, **kwargs)
    
    # Close it
    if close:
        plt.close()

    if verbose:
        print("Done")
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Equations                       ║
# ╚══════════════════════════════════════════════════╝
#region

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
    return psi_max - ((psi_max - psi_min) * term) / (term + 1)

def dX_dt(X, t, psi_max, psi_min, zMIC, k, params, environment):
    '''function in which growth rate is calculated depending on the environmental conditions'''

    # decide in which timestep(e.g day) to quit the administration of antibiotic
    if t > 5: 
        a_t = 3 # antibiotic concentration
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # Environmental variation (as an environmental Cue) at time t
    growth_rate_modifier = psi_max * reaction_norm(params["I0"], params["b"], current_env) # new psimax depending on plasticity
    growth_rate = np.log(10) * psi(a_t, growth_rate_modifier, psi_min, zMIC, k) * X

    return max(growth_rate, -X / 0.04)

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Parameters                      ║
# ╚══════════════════════════════════════════════════╝
#region
# All environments must have different keys otherwise will be overwritten
# All environments must have at least one different value otherwise only the last will be saved
environments_params = {
    "Env 1": {"A": 0.3, "B": 0.0, "L": 10, "R": 2, "t": 110},
    "Env 2": {"A": 0.6, "B": 0.0, "L": 10, "R": 2, "t": 110},
    "Env 3": {"A": 0.9, "B": 0.0, "L": 10, "R": 2, "t": 110},
    "Env 4": {"A": 1.2, "B": 0.0, "L": 10, "R": 2, "t": 110},
    # "Env 5": {"A": 4, "B": 0.0, "L": 10, "R": 2, "t": 110},
}

genotypes_params = {
    "Genotype 1": {"I0": 0.2, "b": 0.8},
    "Genotype 2": {"I0": 0.4, "b":0.6},
    "Genotype 3": {"I0": 0.6, "b": 0.4},
    "Genotype 4": {"I0": 0.8, "b": 0.2},
    "Genotype 5": {"I0": 0.2, "b": 1.4},    
}

psi_min = -2 # maximum death rate
zMIC = 2 # concentration in which net growth rate is zero
k = 0.8  # Using a single mean k value
psi_max = 0.8  # maximal growth rate

time_frame = np.linspace(0, 10, 11) #should be passed on odeint()
initial_populations = [1e3]

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Simulations                     ║
# ╚══════════════════════════════════════════════════╝

# region test simulations

#     #region environment construction
# environment = Environment()
# environment.trim()
# environment.save()
#     #endregion

#     #region norms & responses to environmental variation
# environment.gene_reaction_norms(genotypes_params)
# environment.gene_responses(genotypes_params)
#     #endregion

#     #region bacterial growth simulations
# environment.run_simulation(genotypes_params, initial_populations)
#      #endregion

#endregion

#region main simulations
# simulator = Simulator(environments_params, genotypes_params)
# simulator.run()
#endregion

Determistic = [0.3, 0.6, 0.9]
Stochastic = [0.0,0.1]
Lifespan = [10]
RelativeVariation = [2]
Timesteps = [110]

# number of total envs will be A*B*L*R*t numbers of each parameter so be carefull!
envs = {}
counter = 0
keys =["A", "B", "L", "R", "t"] # list of keys to irerate on dict creation
for A in Determistic:
    for B in Stochastic:
        for L in Lifespan:
            for R in RelativeVariation:
                for t in Timesteps:
                    counter += 1 # counter for dynamic name
                    values = [A, B, L, R, t] # list the parameters to iterate on dict creation
                    params = {keys[i]: values[i] for i in range(len(keys))} # parameters dict creation
                    env = {} # create dict for environment
                    env[f"Env {counter}"] = params # assign parameters dict to env dict
                    print('env:',env)
                    envs.update(env) # add env dict to envrironments dict
print("envs", envs)

print("env params", environments_params)

simulator = Simulator(envs, genotypes_params)
simulator.run()