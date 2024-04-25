# ╔══════════════════════════════════════════════════╗
# ║                   Classes                        ║
# ╚══════════════════════════════════════════════════╝
from tools import *
from equations import *
from parameters import *
from scipy.integrate import odeint


#region

class Environment():

    def __init__(self, A = 0.9, B = 0, L = 10, R = 1, t = 110):
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
        ax.plot(self.t, self.variation, label='Environmental Variation', linestyle="dashdot", color="purple")
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
        None
        # renderPeriod = 3
        # initial_title = self.ax.get_title() # keep original title to put after rendering
        # self.ax.set_title(self.ax.get_title() + "\n" + f'Will disappear after {renderPeriod} seconds') # temp title for rendering
        # plt.show(block=False) # stop from blocking in the execution
        # plt.pause(renderPeriod)
        # self.ax.set_title(initial_title) # set initial title again
        
    def save(self):
        # self.fig.savefig(f"Env. variation t={len(self.t)}, A={self.A}, B={self.B}, L={self.L}, R={self.R}, trimmed={self.trimmed}.png")
        save(f'./report/Env. variation t={len(self.t)}, A={self.A}, B={self.B}, L={self.L}, R={self.R}, trimmed={self.trimmed}.png')

    def gene_responses(self, genotypes):

        fig, ax = self._create_plot() # create a copy of the current variation plot (keep clean the original)

        for name, params in genotypes.items():
            I = reaction_norm(params["I0"], params["b"], self.variation)
            ax.plot(self.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")

        ax.set_title('Phenotypic Responses')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Phenotypic response (I)')
        ax.legend(bbox_to_anchor=(1.34, 1))

        save(f'./report/Responses to Variation')

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

        save(f'./report/Reaction Norms')

    def run_simulation(self, genotypes, initial_populations):

        fig, ax = plt.subplots(figsize=(14,6)) # prepare plot

        # add environmental variation information
        environmental_variation_layer_applier(time_frame,ax,self.variation)

        # add antibiotic exposure information
        antibiotic_exposure_layers_applier(time_frame,ax)

        # plot dynamics
        for initial_population in initial_populations:
            for name, params in genotypes.items():

                X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, self,antibody_concentration)) # args will be passed down to dX_dt
                ax.plot(time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')

        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e9)   

        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.41, 1), fontsize="7") # place legend out of plot

        # fig.savefig(f'Genotypes dynamics.png')
        save(f'./report/Population Dynamics')

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

        self.envs_plot = self.yield_environment_plots_with_antibiotic_frames()
        self.norms_plot = self.yield_reaction_norms()
        self.responses_plot = self.yield_phenotypic_responses()
        self.dynamics_plot = self.yield_population_dynamics_with_antibiotic_frames()

        print("..end of simulation")

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
        # the different values of the parameter of interest become layers (rows/columns) and all the other plots in the layer are scheduled based on 
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
                    row_vector.append(env) # append to the row vector every env that has the A value
            row_vectors.append(row_vector) # append the whole row for that value to the stack

        # based on the number of rows and the length of a row i create dimensions to pass on the add_gridspec()
        rows = len(row_vectors) # turn rows of vectors into rows
        columns = len(row_vectors[0]) # turn length of a row vector into columns

        return row_vectors, rows, columns
    
    def generate_report(self):
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

        fig.add_subplot(221)
        plt.imshow(img1)
        plt.axis('off') 
        fig.add_subplot(222)
        plt.imshow(img2)
        plt.axis('off') 
        fig.add_subplot(223)
        plt.imshow(img3)
        plt.axis('off') 
        fig.add_subplot(224)
        plt.imshow(img4)
        plt.axis('off')

        save('./report/report', dpi=1000)

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
                    custom_plot(axs[row,column],env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}', ylim=(0,1))
                    axs[row,0].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                    axs[-1,column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")

                else:   # if not, it has one dimension
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        custom_plot(axs[row],env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}', ylim=(0,1))
                        axs[row].set_ylabel(f"A ={vector[0].A}", rotation="horizontal", fontsize=14, weight="bold")
                        axs[-1].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")

                    else:               # else another parameter has variation     
                        custom_plot(axs[column],env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}', ylim=(0,1))
                        axs[0].set_ylabel(f"A ={vector[0].A}", rotation="horizontal", fontsize=14, weight="bold")
                        axs[column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")

        
        save('./report/Stacked Environments')
        return fig, axs

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
                        custom_plot(axs[row,column], env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
                else:
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            custom_plot(axs[row], env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                    
                    else:
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            custom_plot(axs[column], env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")

        
        save('./report/Stacked Reaction Norms')
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
                if axs.ndim > 1: # check if grid has two dimensions (+unique values for the another parameter )
                    custom_plot(axs[row,column], env.t, env.variation,  label='Environmental Variation', linestyle="dashdot", color="purple")
                    for name, params in self.genotypes.items():
                        I = reaction_norm(params["I0"], params["b"], env.variation)
                        custom_plot(axs[row,column], env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title=f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
                else:
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        custom_plot(axs[row], env.t, env.variation,  label='Environmental Variation', linestyle="dashdot", color="purple")
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            custom_plot(axs[row], env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title=f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                           
                    else:
                        custom_plot(axs[column], env.t, env.variation,  label='Environmental Variation', linestyle="dashdot", color="purple")
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            custom_plot(axs[column], env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title=f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                                                      

        save('./report/Stacked Phenotypics Responses')
        return fig
    
    def yield_population_dynamics(self):

        row_vectors, rows, columns = self._plot_layer_constructor()

        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(("Population Dynamics \n"
                        r"$\bf{" + "Response Curve Parameters:" + "}$"
                        f"k={k}, "
                        f"Ψmax={psi_max}, "
                        f"Ψmin={psi_min}, "
                        f"MIC={zMIC}, "
                        f"c={antibody_concentration} \n"
                        r"$\bf{" + "Growth  Rate  Modifier:" + "}$" + f"{get_function_body(growth_rate_modifier)} \n"
                        r"$\bf{" + "Death  Rate  Modifier:" + "}$" + f"{get_function_body(death_rate_modifier)}"
                        ),
                        fontsize=20)
        fig.text(0.5, 0.07, "Time (t)", fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Bacterial density", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        for row, vector in enumerate(row_vectors):
            for column, env in enumerate(vector):
                if axs.ndim > 1: # check if grid has two dimensions (+unique values for the another parameter )
                    
                    for initial_population in initial_populations:
                        for name, params in self.genotypes.items():
                            X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibody_concentration)) # args will be passed down to dX_dt
                            custom_plot(axs[row,column], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))
                            axs[row,0].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                            axs[-1,column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")                       
                else:
                    if len(row_vectors)>1: # check if the the parameter of interest has more than one value
                        for initial_population in initial_populations:
                            for name, params in self.genotypes.items():
                                X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibody_concentration)) # args will be passed down to dX_dt
                                custom_plot(axs[row], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))
                                axs[row].set_ylabel(f"A ={vector[0].A}", rotation="horizontal", fontsize=14, weight="bold")
                                axs[-1].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")                
                    
                    else:   # else another parameter has variation 
                         for initial_population in initial_populations:
                            for name, params in self.genotypes.items():
                                X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibody_concentration)) # args will be passed down to dX_dt
                                custom_plot(axs[column], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))                                
                                axs[0].set_ylabel(f"A ={vector[0].A}", rotation="horizontal", fontsize=14, weight="bold")
                                axs[column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")
        
        save('./report/Stacked Population Dynamics')
        return fig, axs
    
    def yield_environment_plots_with_antibiotic_frames(self):
        # this function works right but has problem with the custom save function & that is why i dont use it.
        # problem is that custom save() doesnt take this fig as current and it saves an irrelevant plot of the past
        # when tried to activate fig with plt.figure(fig) saving it causes a traceback error with save() (somethings broken with plt.close())
        

        fig, axs = self.yield_environment_plots()
        
        if axs.ndim > 1:
            for vector in axs:
                for ax in vector:
                    antibiotic_exposure_layers_applier(time_frame,ax)
        else:
            for ax in axs:
                antibiotic_exposure_layers_applier(time_frame,ax)

        fig.savefig("./report/Stack environments with Antibiotics Layers")
        return fig
        # fig.savefig('./report/stacked environments with Antibiotc layers')
        # pl.figure(figure)
        # # plt.close()
        # save('./report/Stacked Environments with Antibiotic Layers')
        # print("after")
        
    def yield_population_dynamics_with_antibiotic_frames(self):
        # this function works right but has problem with the custom save function & that is why i dont use it.
        # problem is that custom save() doesnt take this fig as current and it saves an irrelevant plot of the past
        # when tried to activate fig with plt.figure(fig) saving it causes a traceback error with save() (somethings broken with plt.close())

        fig, axs = self.yield_population_dynamics()

        if axs.ndim > 1:
            for vector in axs:
                for ax in vector:
                    antibiotic_exposure_layers_applier(time_frame,ax)
        else:
            for ax in axs:
                antibiotic_exposure_layers_applier(time_frame,ax)

        fig.savefig("./report/Stack Population Dynamics with Antibiotics Layers")
        return fig, axs

    def yield_population_dynamics_with_antibiotic_frames_env_variation(self):
        
        row_vectors, rows, columns = self._plot_layer_constructor()

        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(("Population Dynamics \n"
                        r"$\bf{" + "Response Curve Parameters:" + "}$"
                        f"k={k}, "
                        f"Ψmax={psi_max}, "
                        f"Ψmin={psi_min}, "
                        f"MIC={zMIC}, "
                        f"c={antibody_concentration} \n"
                        r"$\bf{" + "Growth  Rate  Modifier:" + "}$" + f"{get_function_body(growth_rate_modifier)} \n"
                        r"$\bf{" + "Death  Rate  Modifier:" + "}$" + f"{get_function_body(death_rate_modifier)}"
                        ),
                        fontsize=20)
        fig.text(0.5, 0.07, "Time (t)", fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Bacterial density", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        for row, vector in enumerate(row_vectors):
            for column, env in enumerate(vector):
                if axs.ndim > 1: # check if grid has two dimensions (+unique values for the another parameter )
                    
                    for initial_population in initial_populations:
                        for name, params in self.genotypes.items():
                            X = odeint(dX_dt, initial_population, time_frame, 
                                       args=(psi_max, psi_min, zMIC, k, params, env, antibody_concentration)) # args will be passed down to dX_dt
                            custom_plot(axs[row,column], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))
                            axs[row,0].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                            axs[-1,column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")                       
                    
                    # add environmental variation information
                    environmental_variation_layer_applier(time_frame,axs[row,column],env.variation)
                    
                    # add antibiotic exposure information
                    antibiotic_exposure_layers_applier(time_frame,axs[row,column])

                else:
                    if len(row_vectors)>1: # check if the the only parameter of interest has more than one value

                        # add environmental variation information
                        environmental_variation_layer_applier(time_frame,axs[row],env.variation)
                        
                        # add antibiotic exposure information
                        antibiotic_exposure_layers_applier(time_frame,axs[row])
                        
                        for initial_population in initial_populations:
                            for name, params in self.genotypes.items():
                                X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibody_concentration)) # args will be passed down to dX_dt
                                custom_plot(axs[row], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))
                                axs[row].set_ylabel(f"A ={vector[0].A}", rotation="horizontal", fontsize=14, weight="bold")
                                axs[-1].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")                
                    
                    else:   # else another parameter has variation 
                        
                        # add environmental variation information
                        environmental_variation_layer_applier(time_frame,axs[column],env.variation)
                                         
                        # add antibiotic exposure information
                        antibiotic_exposure_layers_applier(time_frame,axs[column]) 
                         
                        for initial_population in initial_populations:
                            for name, params in self.genotypes.items():
                                X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibody_concentration)) # args will be passed down to dX_dt
                                custom_plot(axs[column], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))                                
                                axs[0].set_ylabel(f"A ={vector[0].A}", rotation="horizontal", fontsize=14, weight="bold")
                                axs[column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")

        save('./report/Population Dynamics & antibiotic frames & env variation')
        return fig, axs
#endregion