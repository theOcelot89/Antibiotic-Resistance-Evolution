# ╔══════════════════════════════════════════════════╗
# ║                   Imports                        ║
# ╚══════════════════════════════════════════════════╝
#region
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from scipy.integrate import odeint
import io 
import os
from PIL import Image 
import pylab as pl
import inspect
import ast
from parameters import *
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                   Classes                        ║
# ╚══════════════════════════════════════════════════╝
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

def antibiotic_exposure_layers_applier(period, ax):

    #create vectors for the different situations you wish to highlight
    antibiotic_exposure_frame = [] 
    antibiotic_NOT_exposure_frame =[]

    for time in period:
        if is_time_for_administration(time) :
            antibiotic_exposure_frame.append(int(time)) 
        else:
            antibiotic_NOT_exposure_frame.append(int(time))

    # appending highlight to the plot
    for i in set(antibiotic_exposure_frame):
        ax.axvspan(i, i+1, facecolor='lightcoral', edgecolor='none', alpha=0.3 ) 
    for i in set(antibiotic_NOT_exposure_frame):
        ax.axvspan(i, i+1, facecolor='palegreen', edgecolor='none', alpha=0.3 )

    #create color patches for the legend to show
    exposure_patch = mpatches.Patch(color='red',  alpha=.2, label='Antibiotic Exposure')
    no_exposure_patch = mpatches.Patch(color='green', alpha=.2, label='No exposure')

    # https://www.statology.org/matplotlib-manual-legend/
    handles, labels = ax.get_legend_handles_labels() # extracting the previous legend stuff
    handles.extend([exposure_patch,no_exposure_patch]) # adding the patch to the old stuff
    ax.legend(handles=handles)

    return ax

def environmental_variation_layer_applier(time_frame, ax, variation):

    if len(variation) < int(max(time_frame))+1:
        raise Exception("your time frame is is bigger than time. Please reduce time frame or increase time..")
    
    variation_axe = ax.twinx()

    # here we use max value of time frame for time reference
    # but we add 1 in order to index correctly the variation list until the right time 
    # and because x and y must be of the same length we have to raise also the the_frame reference  
    custom_plot(variation_axe, np.arange(int(max(time_frame))+1), variation[:int(max(time_frame))+1], linestyle="dashdot", color="purple", alpha=0.3, ylim=(0,1))
    variation_axe.yaxis.set_major_locator(ticker.NullLocator()) # remove ticks and labels rom y axis




def custom_plot(ax, xdim, ydim, **params):

    # extracting necessary parameters for plot() method or set default value to None or others
    linestyle = params.get('linestyle', None) 
    color = params.get('color', None) 
    label = params.get('label', None)
    alpha = params.get('alpha', None)

    ax.plot(xdim, ydim, linestyle=linestyle, color=color, label=label, alpha=alpha)

    if "legend" in params:
        ax.legend()

    if "legend_title" in params:
        ax.legend(title=params["legend_title"])

    if "ylim" in params:
        ax.set_ylim(params["ylim"][0],params["ylim"][1])

    if "yscale" in params:
        ax.set_yscale(params["yscale"])

def get_function_body(func):
    # Get the source code of the function as a string
    source_code = inspect.getsource(func)
    
    # Parse the source code into an Abstract Syntax Tree (AST)
    tree = ast.parse(source_code)
    
    # Extract the body of the function
    function_body = tree.body[0].body
    
    # Convert the body back into a string
    function_body_str = ast.unparse(function_body)

    # Remove the "return" statement if present
    if function_body_str.strip().startswith("return "):
        function_body_str = function_body_str.replace("return ", "", 1)
    
    return function_body_str

def bold_text(text):
  return "\033[1m" + text + "\033[0m"
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
    return (psi_max - psi_min) * (term / (term - psi_min/psi_max))
    return psi_max - ((psi_max - psi_min) * term) / (term + 1) # Giorgio's Implementation

def dX_dt(X, t, psi_max, psi_min, zMIC, k, params, environment,antibody_concentration):
    '''function in which growth rate is calculated depending on the environmental conditions'''

    if population_is_below_threshold(X,10):
        X = 0

    if is_time_for_administration(t): 
        a_t = antibody_concentration 
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # Environmental variation (as an environmental Cue) at time t
    modified_current_env = current_env * (1 - (X/1e9))

    modified_growth_rate = growth_rate_modifier(psi_max, params, modified_current_env)
    modified_death_rate = death_rate_modifier(modified_growth_rate)
    growth_rate_after_antibiotic = modified_growth_rate -  psi(a_t, modified_growth_rate, modified_death_rate, zMIC, k)
    actual_growth_rate = np.log(10) * growth_rate_after_antibiotic * X * (1 - (X/1e9))

    return max(actual_growth_rate, -X / 0.04)

def is_time_for_administration(time):
    # not statement reverses the antibiotic exposure time frames (simply put in front of expression)
    return time % 20 < 10

# def is_time_for_delution(time):
#     return time % 10 < 3

def population_is_below_threshold(X, threshold):
    return X < threshold


def growth_rate_modifier(psi_max, params, env):
    return psi_max * reaction_norm(params["I0"], params["b"], env)

def death_rate_modifier(growth):
    return  - growth * 1.5
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Parameters                      ║
# ╚══════════════════════════════════════════════════╝
#region
# All environments must have different keys otherwise will be overwritten
# All environments must have at least one different value otherwise only the last will be saved
# determistic = [0.3,0.6]
# stochastic = [0.0,]
# lifespan = [10]
# relativeVariation = [1,2]
# timesteps = [101]

# environments_params = construct_params(determistic, stochastic, lifespan, relativeVariation, timesteps)


# genotypes_params = {
#     "Genotype 1": {"I0": 0.1, "b": 0.9},
#     # "Genotype 2": {"I0": 0.4, "b":0.6},
#     "Genotype 3": {"I0": 0.5, "b": 0.4},
#     "Genotype 4": {"I0": 0.8, "b": 0},
#     # "Genotype 5": {"I0": 0.2, "b": 1.4},    
# }

# antibody_concentration = 100
# psi_min = -2 # maximum death rate
# zMIC = 2 # concentration in which net growth rate is zero
# k = 0.8  # Using a single mean k value
# psi_max = 0.3  # maximal growth rate
# initial_populations = [1e7]



# # for all simulations and layer appliers to work properly
# # the slicing must be at least time+1 (e.g. 101 slices for time=100)
# time_frame = np.linspace(0, 100, 101) #should be passed on odeint()



#endregion

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







