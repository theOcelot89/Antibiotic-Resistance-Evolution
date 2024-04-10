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
        
        
        print(f"New environment created!, params: A={self.A}, B={self.B}, L={self.L}, R={self.R}, t={len(self.t)} ")

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
        print(f"simulator contains {len(self.environments)} environments & {len(self.genotypes)} genotypes.")

         # code will not work with one environment, so i put this condition
        if len(self.environments) == 0 or len(self.environments) == 1:
            print("please provide more than 1 environment")
            exit()
    
    def run(self):

        print("Simulation begins...")

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

    def _plot_layer_constructor(self):
        # here i make a method of constructing layers of plots based on a parameter of interest (here is A)
        # the different values of the parameter of interest become layers and all the other plots in the layer are scheduled based on 
        # on these layers (hierarchy on the other axis is based on the line of the environments that are feeded and this 
        # functionality is created internally - i didnt code for this)
        # IMPORTANT!!! rows here is just a name convension. The code for plotting the layers horizontally/vertically
        # is based on how you pass row/columns in add_gridspec() and axs[i,j] OR axs[j,i]
        # you have to reverse these parameters in order to change direction of plotting.
        A_list  = []  
        for env in self.environments:  
            A_list.append(env.A)

        # here i construct row vectors (which means i ll plot in row layers)
        row_vectors = [] # stack of rows
        for A in sorted(set(A_list)): # sorted set of unique values of A
            row_vector = [] # construct a row vector for every unique A variable 
            for env in self.environments:
                if A == env.A:
                    row_vector.append(env) # append to the row vector every env that has the A parameter
            row_vectors.append(row_vector)

        # based on the number of rows and the length of a row i create dimensions to pass on the add_gridspec()
        rows = len(row_vectors) # turn rows of vectors into rows
        columns = len(row_vectors[0]) # turn length of a row vector into columns

        return row_vectors, rows, columns
    
    def _generate_report(self):
        # ok that was a touch one & a number of steps & guides have to be considered in order for a nice report sheet
        # 1. convert figures to images (that is the only way to combine them) https://www.geeksforgeeks.org/saving-a-plot-as-an-image-in-python/
        # 2. use fig.figAddSubplot technique to create the grid https://www.geeksforgeeks.org/how-to-display-multiple-images-in-one-figure-correctly-in-matplotlib/
        # 3. use plt.imshow to render the image in the figure https://www.geeksforgeeks.org/how-to-display-multiple-images-in-one-figure-correctly-in-matplotlib/
         
        # print(self.norms_plot.get_size_inches())  #dynamically set width/height from the dimensions of the plots

        img1 = fig2img(self.dynamics_plot)
        img2 = fig2img(self.responses_plot)
        img3 = fig2img(self.norms_plot)
        img4 = fig2img(self.envs_plot)

        fig = plt.figure(layout="compressed")

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

        # fig.savefig('Report', dpi=600, bbox_inches='tight') # dpi for a better resolution, bboxinches for trimming margins
        save('./report/report', dpi=6000)

    def yield_environment_plots(self):
        # this technique is based on matplots basic tutorial
        # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
        
        row_vectors, rows, columns = self._plot_layer_constructor()

        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(u"Environmental Variations\nE\u209C = A·sin(2πt/LR) + B·ε", fontsize = 30, y= 0.95)
        fig.text(0.5, 0.07, "Time (t)", va="center",  fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Environmental variation (E)", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        for row, vector in enumerate(row_vectors): # iterate on the stack of rows
            for column, env in enumerate(vector): # iterate on the row itself
                if axs.ndim > 1: # check if grid has two dimensions (+unique values for the another parameter )
                    axs[row,column].plot(env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}')
                    axs[row,column].legend()
                    axs[row,column].set_ylim(0,1)
                else:   # if not, it has one dimension
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        axs[row].plot(env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}')
                        axs[row].legend()
                        axs[row].set_ylim(0,1)
                    else:               # else another parameter has variation     
                        axs[column].plot(env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}')
                        axs[column].legend()
                        axs[column].set_ylim(0,1)

        save('./report/Stacked Environments')
        print("environment plots DONE")
        return fig

    def yield_reaction_norms(self):

        row_vectors, rows, columns = self._plot_layer_constructor()       
        
        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0.1, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.text(0.35, 0.07, "Environmental variation (E)",   fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Response (I)", rotation="vertical", va="center", fontsize=20) # put only 1 y label
        fig.suptitle(u"Reaction Norms\n I=I\u2080 + b·C", fontsize = 30, y= 0.95)
        
        for row, vector in enumerate(row_vectors):
            for column, env in enumerate(vector):
                if axs.ndim > 1: # check if grid has two dimensions (+unique values for the another parameter )
                    for name, params in self.genotypes.items():
                        I = reaction_norm(params["I0"], params["b"], env.variation)
                        axs[row,column].plot(env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                        axs[row,column].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
                else:
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            axs[row].plot(env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                            axs[row].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                            
                    else:
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            axs[column].plot(env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                            axs[column].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                                                        


        save('./report/Stacked Reaction Norms')
        print("reaction norms DONE")
        return fig
    
    def yield_phenotypic_responses(self):
        
        row_vectors, rows, columns = self._plot_layer_constructor()

        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(u"Phenotypic Responses\n I\u209C=I\u2080 + b·E\u209C", fontsize = 30, y= 0.95)
        fig.text(0.5, 0.07, "Time (t)",   fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Response (I)", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        for row, vector in enumerate(row_vectors):
            for column, env in enumerate(vector):
                axs[row,column].plot(env.t, env.variation, label='Environmental Variation', linestyle="dashdot")
                if axs.ndim > 1: # check if grid has two dimensions (+unique values for the another parameter )
                    for name, params in self.genotypes.items():
                        I = reaction_norm(params["I0"], params["b"], env.variation)
                        axs[row,column].plot(env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                        axs[row,column].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
                else:
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            axs[row].plot(env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                            axs[row].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                            
                    else:
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            axs[column].plot(env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")
                            axs[column].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                                                        

        save('./report/Stacked Phenotypics Responses')
        print("phenotypics responses DONE")
        return fig
    
    def yield_population_dynamics(self):

        fig = plt.figure(figsize=(12, len(self.environments)*6)) # empty figure for template, dynamic height of plot
        gs = fig.add_gridspec(len(self.environments), hspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(f"Population Dynamics \nResponse Curve Parameters: k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}", fontsize = 20, y= 0.95)
        fig.text(0.5, 0.07, "Time (t)", fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Bacterial density", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        for i, env in enumerate(self.environments):
            for initial_population in initial_populations:
                for name, params in self.genotypes.items():

                    X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env)) # args will be passed down to dX_dt
                    axs[i].plot(time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}')
                    axs[i].set_yscale('log')
                    axs[i].set_ylim(1, 1e9) 
                    axs[i].legend(title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")  

        save('./report/Stacked Population Dynamics')
        print("population dynamics DONE")
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

def construct_params(determistic, stochastic, lifespan, relativeVariation, timesteps):

    # number of total envs will be A*B*L*R*t of different parameter so be carefull!
    envs = {}
    counter = 0
    keys =["A", "B", "L", "R", "t"] # list of keys to irerate on dict creation

    for A in determistic:
        for B in stochastic:
            for L in lifespan:
                for R in relativeVariation:
                    for t in timesteps:
                        counter += 1 # counter for dynamic env name
                        values = [A, B, L, R, t] # list the parameters to iterate on dict creation
                        params = {keys[i]: values[i] for i in range(len(keys))} # parameters dict creation
                        env = {} # create dict for environment
                        env[f"Env {counter}"] = params # assign parameters dict to env dict
                        envs.update(env) # add env dict to envrironments dict
    
    return envs
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
    # "Env 1": {"A": 0.3, "B": 0.0, "L": 10, "R": 2, "t": 110},
    # "Env 2": {"A": 0.6, "B": 0.0, "L": 10, "R": 2, "t": 110},
    # "Env 3": {"A": 0.9, "B": 0.0, "L": 10, "R": 2, "t": 110},
    # "Env 4": {"A": 1.2, "B": 0.0, "L": 10, "R": 2, "t": 110},
    # "Env 5": {"A": 4, "B": 0.0, "L": 10, "R": 2, "t": 110},
}

determistic = [0.3, 0.6, 0.9, 1]
stochastic = [0.0]
lifespan = [10]
relativeVariation = [2,4,6]
timesteps = [110]

environments_params = construct_params(determistic, stochastic, lifespan, relativeVariation, timesteps)


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
simulator = Simulator(environments_params, genotypes_params)
# simulator.run()
simulator.yield_environment_plots()
simulator.yield_phenotypic_responses()
simulator.yield_reaction_norms()
#endregion






