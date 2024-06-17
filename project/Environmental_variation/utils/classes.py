# ╔══════════════════════════════════════════════════╗
# ║                   Classes                        ║
# ╚══════════════════════════════════════════════════╝
from matplotlib.lines import Line2D

from .equations import *
from .tools import *
from scipy.integrate import odeint, solve_ivp   


if __name__ == "__main__":
    print('classes called directly nothing to show..')

else:
    print('classes loaded..')

    
class Environment():

    def __init__(self, A , B , L , R , t, genotypes, framework):
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
        self.genotypes = genotypes
        self.framework = framework 

        self.env_params = self.A, self.B, self.L, self.R
        self.results, self.fig, self.ax = self._simulation_odeint_with_mutation()

    def _simulation(self):

        env_params = self.env_params
        genotypes = self.genotypes
        framework = self.framework

        initial_populations = framework["Initial Populations"]
        time_frame = framework["time frame"]        
        
        results = {}
        for initial_population in initial_populations:
            y0 = [initial_population,0,0,0,0,0,0,0]
            for name, params in genotypes.items():
                X = odeint(sim, y0, time_frame,args=(env_params, params, framework)) 
                results[name] = X
        return results

    def _simulation_odeint_with_mutation(self):

        env_params = self.env_params
        genotypes = self.genotypes
        framework = self.framework

        initial_populations = framework["Initial Populations"]
        mutation_population = 0 # placeholder for mutants to induce into the simulation
        time_frame = framework["time frame"]         

        results = {}
        for initial_population in initial_populations:
            for name, params in genotypes.items():
                
                y0 = [initial_population, mutation_population, 0, 0, 0, 0, 0, 0, 0]
                X = odeint(sim_with_mutation_event, y0, time_frame, args=(env_params, params, framework)) 
                results[name] = X


        # PLOT THE RESULTS
        fig , ax = plt.subplots(figsize=(18,6))

        for name, result in results.items():

            wild_type_dynamics = result[:,0]
            mutant_type_dynamics = result[:,1]
            
            ax.plot(framework['time frame'], wild_type_dynamics, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")
            ax.plot(framework['time frame'], mutant_type_dynamics, label=f"mutant: {name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")


        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e10)                   
        ax.legend()

        # PLACE LEGEND OUT OF PLOT
        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot 


        save("./results/Dynamics with mutation", close=False)
        return results, fig, ax

    def _simulation_with_mutation_event(self):

        env_params = self.env_params
        genotypes = self.genotypes
        framework = self.framework

        initial_populations = framework["Initial Populations"]
        mutation_population = 0 # placeholder for mutants to induce into the simulation
        time_frame = framework["time frame"]   
        start = (int(time_frame[0])) # index time frame for solve_ivp
        end = (int(time_frame[-1]))  # index time frame for solve_ivp
        t_span = np.array([start,end]) # time span that is mandatory for solve_ivp  

        results = {}
        for initial_population in initial_populations:
            for name, params in genotypes.items():
                
                y0 = np.array([initial_population, mutation_population, 0, 0, 0, 0, 0, 0, 0])
                X = solve_ivp(sim_mutation_VERSION2, t_span, y0,  args=(env_params, params, framework), t_eval= time_frame, method='LSODA') 
                results[name] = X


        # PLOT THE RESULTS
        fig , ax = plt.subplots(figsize=(14,6))

        for name, result in results.items():

            time = result.t
            wild_type_dynamics = result.y[0]
            mutant_type_dynamics = result.y[1]
            
            ax.plot(time, wild_type_dynamics, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")
            ax.plot(time, mutant_type_dynamics, label=f"mutant: {name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")


        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e10)                   
        ax.legend()
        save("./results/Dynamics with mutation VERSION 2", close=False)
        return fig, ax

    def _simulation_with_event(self):


        env_params = self.env_params
        genotypes = self.genotypes
        framework = self.framework

        initial_populations = framework["Initial Populations"]
        mutation_population = 0 # placeholder for mutants to induce into the simulation
        time_frame = framework["time frame"]   
        start = (int(time_frame[0])) # index time frame for solve_ivp
        end = (int(time_frame[-1]))  # index time frame for solve_ivp
        t_span = np.array([start,end]) # time span that is mandatory for solve_ivp  

                
        def event(t,y,env_params, params, framework):
            return int(t) - 150 # the return needs to be 0 for the event to trigger
        # event.terminal = True

        results = {}
        for initial_population in initial_populations:
            for name, params in genotypes.items():
                
                y0 = np.array([initial_population, 0, 0, 0, 0, 0, 0, 0])
                X = solve_ivp(sim_ivp, t_span, y0,  
                              args=(env_params, params, framework), 
                              t_eval= time_frame, 
                              method='LSODA',
                              events=event,
                              dense_output=True) 
                results[name] = X

                # print(X.t_events)
                # print(X.sol(X.t_events[0][0]))
                # print("x.y:", X.y)
        

        # PLOT THE RESULTS
        fig , ax = plt.subplots(figsize=(14,6))

        for name, result in results.items():

            time = result.t
            wild_type_dynamics = result.y[0]
            # mutant_type_dynamics = result.y[1]
            print('time events', result.t_events)
            print('y events', result.y_events)
            
            ax.plot(time, wild_type_dynamics, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")
            ax.plot(result.t_events[0][0], result.y_events[0][0][0], 'o')
            # ax.plot(time, mutant_type_dynamics, label=f"mutant: {name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")


        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e10)                   
        ax.legend()
        save("./results/Dynamics with event", close=False)
        return fig, ax

    def variation(self):

        time_frame = self.framework["time frame"]

        index = list(self.results)[0] # grab the name of the first key in result in order to index in the next step
        X = self.results[index] # use the index to take the first genotype results (here we dont care about which genotypes as "true variation" is independent from it.)
        params = self.genotypes[index] # also use the index to grab the params of the genotype (also here true variation is not affected from this)
        env_params = self.env_params
        framework = self.framework

        # https://stackoverflow.com/questions/54365358/odeint-function-from-scipy-integrate-gives-wrong-result
        # i use this code in order to draw the true variation information that i want in order to plot correctly
        variation = [sim(y, time, env_params, params, framework)[1] for time, y in zip(time_frame, X)]

        fig , ax = plt.subplots(figsize=(14,6))
        ax.plot(time_frame, variation, linestyle= "dashdot", color="purple", label="True Variation")
        ax.legend()
        ax.grid()

        # # PLACE LEGEND OUT OF PLOT
        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot

        self.variation = variation
        save('./results/Environmental Variation')

    def normalized_variation(self):

        time_frame = self.framework["time frame"]

        index = list(self.results)[0] # grab the name of the first key in result in order to index in the next step
        X = self.results[index] # use the index to take the first genotype results (here we dont care about which genotypes as "true variation" is independent from it.)
        params = self.genotypes[index] # also use the index to grab the params of the genotype (also here true variation is not affected from this)
        env_params = self.env_params
        framework = self.framework

        # https://stackoverflow.com/questions/54365358/odeint-function-from-scipy-integrate-gives-wrong-result
        # i use this code in order to draw the true variation information that i want in order to plot correctly
        variation = [sim(y, time, env_params, params, framework)[7] for time, y in zip(time_frame, X)]

        fig , ax = plt.subplots(figsize=(14,6))
        ax.plot(time_frame, variation, linestyle= "dashdot", color="purple", label="True Variation")
        ax.legend()
        ax.grid()

        self.variation = variation
        save('./results/Normalized Env Variation')

    def responses(self):

        time_frame = self.framework["time frame"]

        results = self.results
        genotypes = self.genotypes
        framework = self.framework
        env_params = self.env_params
        fig , ax = plt.subplots(figsize=(14,6))

        for name, X in results.items():
            # https://stackoverflow.com/questions/54365358/odeint-function-from-scipy-integrate-gives-wrong-result
            # i use this code in order to draw the true variation information that i want in order to plot correctly
            response = [sim(y, time, env_params, genotypes[name], framework)[2] for time, y in zip(time_frame, X)]
            ax.plot(time_frame, response, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")

        ax.plot(time_frame, self.variation, linestyle= "dashdot", color="purple", label="True Variation")
        ax.set_title('Phenotypic Responses')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Phenotypic response (I)')          
        ax.legend()

        # # PLACE LEGEND OUT OF PLOT
        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot

        save("./results/Responses")

    def actual_psi_max(self):

        time_frame = self.framework["time frame"]

        results = self.results
        genotypes = self.genotypes
        framework = self.framework
        env_params = self.env_params
        color_list = generate_color_list(len(genotypes)) 

        fig , ax = plt.subplots(figsize=(14,6))

        for index,(name, X) in enumerate(results.items()):
            # https://stackoverflow.com/questions/54365358/odeint-function-from-scipy-integrate-gives-wrong-result
            # i use this code in order to draw the true variation information that i want in order to plot correctly

            # draw psi max/min data from the results and plot
            psi_max = [sim(y, time, env_params, genotypes[name], framework)[3] for time, y in zip(time_frame, X)]
            psi_min = [sim(y, time, env_params, genotypes[name], framework)[4] for time, y in zip(time_frame, X)]
            ax.plot(time_frame, psi_max, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}", color=color_list[index])
            ax.plot(time_frame, psi_min, label=f" psi min  {name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}", color=color_list[index], linestyle="dashdot")

        ax.set_title(f'Actual Ψmax = {get_function_body(growth_rate_modifier)}, Ψmin = {get_function_body(death_rate_modifier)}')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Ψmax')          
        ax.legend()
        save("./results/Actual psi Max")
    
    def actual_psi_max__antibiotic_effect__growth_rate(self):

        time_frame = self.framework["time frame"]

        results = self.results
        genotypes = self.genotypes
        framework = self.framework
        env_params = self.env_params

        for name, X in results.items():

            fig , ax = plt.subplots(figsize=(14,6))

            # https://stackoverflow.com/questions/54365358/odeint-function-from-scipy-integrate-gives-wrong-result
            # i use this code in order to draw the true variation information that i want in order to plot correctly
            modified_psi_max = [sim(y, time, env_params, genotypes[name], framework)[3] for time, y in zip(time_frame, X)]
            antibiotic_effect = [sim(y, time, env_params, genotypes[name], framework)[6] for time, y in zip(time_frame, X)]
            growth = [sim(y, time, env_params, genotypes[name], framework)[5] for time, y in zip(time_frame, X)]
            ax.plot(time_frame, modified_psi_max, label=f"Modified Ψmax", linestyle="dashdot")
            ax.plot(time_frame, antibiotic_effect, label=f"Antibiotic effect", linestyle="dashdot")
            ax.plot(time_frame, growth, label=f"Growth Rate (Ψmax -  Antibiotic effect)")           
            ax.set_title(f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")
            ax.set_xlabel('Time (t)')
            ax.set_ylabel('Rate')  
            ax.set_ylim(-1,1)        
            ax.legend()

            antibiotic_exposure_layers_applier(time_frame, ax)

            # PLACE LEGEND OUT OF PLOT
            pos = ax.get_position() #returns bbox in order to manipulate width/height
            ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
            ax.legend(bbox_to_anchor=(1.02, 1)) # place legend out of plot

            save(f"./results/{name}, I0 {genotypes[name]["I0"]}, b {genotypes[name]["b"]}")

    def growth_rate_after_antibiotic(self):

        time_frame = self.framework["time frame"]

        results = self.results
        genotypes = self.genotypes
        framework = self.framework
        env_params = self.env_params
        fig , ax = plt.subplots(figsize=(14,6))

        for name, X in results.items():
            # https://stackoverflow.com/questions/54365358/odeint-function-from-scipy-integrate-gives-wrong-result
            # i use this code in order to draw the true variation information that i want in order to plot correctly
            response = [sim(y, time, env_params, genotypes[name], framework)[5] for time, y in zip(time_frame, X)]
            ax.plot(time_frame, response, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")

        ax.set_title('Growth Rate after antibiotic effect')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel(' Growth rate')          
        ax.legend()

        antibiotic_exposure_layers_applier(time_frame, ax)
        save("./results/Growth Rate")

    def dynamics(self):

        results = self.results
        genotypes = self.genotypes
        framework = self.framework

        fig , ax = plt.subplots(figsize=(14,6))

        for name, result in results.items():

            dynamics = result[:,0]
            ax.plot(framework['time frame'], dynamics, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")

        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e10)                   
        ax.legend()
        save("./results/Dynamics", close=False)
        return fig, ax

    def dynamics_with_mutation(self):

            results = self.results
            genotypes = self.genotypes
            framework = self.framework

            fig , ax = plt.subplots(figsize=(14,6))

            for name, result in results.items():

                wild_type_dynamics = result[:,0]
                mutant_type_dynamics = result[:,1]
                
                ax.plot(framework['time frame'], wild_type_dynamics, label=f"{name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")
                ax.plot(framework['time frame'], mutant_type_dynamics, label=f"mutant: {name}, I0:{genotypes[name]["I0"]}, b:{genotypes[name]["b"]}")


            ax.set_xlabel('Time')
            ax.set_ylabel('Bacterial Density')
            ax.set_yscale('log')
            ax.set_ylim(1, 1e10)                   
            ax.legend()
            save("./results/Dynamics with mutation", close=False)
            return fig, ax
    
    def dynamics_with_antibiotic_frames(self):

        time_frame = self.framework["time frame"]        
        fig, ax = self.fig, self.ax

        # add antibiotic exposure information
        antibiotic_exposure_layers_applier(time_frame,ax)

        # # PLACE LEGEND OUT OF PLOT
        # pos = ax.get_position() #returns bbox in order to manipulate width/height
        # ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot 

        save('./results/Dynamics & Antibiotic Frames', close=False)

        self.fig, self.ax = fig, ax
        return fig, ax

    def dynamics_with_antibiotic_frames_and_variation(self):
        
        time_frame = self.framework["time frame"]        
        variation = self.variation
        fig, ax = self.fig, self.ax

        # # PLACE LEGEND OUT OF PLOT
        # pos = ax.get_position() #returns bbox in order to manipulate width/height
        # ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        # ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot
        
        # add environmental variation information
        environmental_variation_layer_applier(time_frame, ax, variation)

        save('./results/Dynamics & Antibiotic Frames & Variation', close=False)
        return fig, ax

    def _create_plot(self):
            
        fig, ax = plt.subplots(figsize=(12, 6))   
        ax.plot(self.t, self.variation, label='Environmental Variation', linestyle="dashdot", color="purple")
        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 
        ax.set_title(f' Environmental variation A={self.A}, B={self.B}, L={self.L}, R={self.R}, trimmed={self.trimmed} ')
        ax.grid()
        return fig, ax
        
    def trim(self):
        '''limiting environmental variation between 0 & 1'''
        # because trim will work directly on self.variation its action is irreversible & will follow the instance for the rest of the script when plots are constructed 

        self.variation[self.variation < 0] = 0  # if any of the values are negative, set them to 0
        self.variation[self.variation > 1] = 1    # if any of the values are greater than 1, set them to 1
        self.trimmed = True

        # reconstruct variation plot
        self.fig, self.ax = self._create_plot()

    def gene_responses(self, genotypes):

        fig, ax = self._create_plot() # create a copy of the current variation plot (keep clean the original)

        for name, params in genotypes.items():
            I = reaction_norm(params["I0"], params["b"], self.variation)
            ax.plot(self.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}")

        ax.set_title('Phenotypic Responses')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Phenotypic response (I)')
        ax.legend(bbox_to_anchor=(1.34, 1))

        save(f'./report/Responses to Variation', close=False)

        return fig

    def gene_reaction_norms(self):

        genotypes = self.genotypes
        variation = self.variation
        variation = np.array(variation) # convert list to numpy array for the reaction norm function
        fig, ax = plt.subplots(figsize=(12,6)) # create plot from scratch

        for name, params in genotypes.items():
            I = reaction_norm(params, variation)
            ax.plot(variation, I, label=f"{name}, I0={params["I0"]}, b={params["b"]}")

        pos = ax.get_position() #returns bbox in order to manipulate width/height
        ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        ax.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot 

        ax.set_title('Reaction Norms')
        ax.set_xlabel('Environmental Variation (E)')
        ax.set_ylabel('Phenotypic response (I)')

        # if not is_called_from_another_function():
        save(f'./results/Norms', close = False)

        return fig

   

        # Unpacking the dictionary into variables
        zMIC, antibiotic_concentration, psi_max, psi_min, k, time_frame, initial_populations = (
        antibiotic_framework["zMIC"],
        antibiotic_framework["Antibiotic Concentration"],
        antibiotic_framework["psi_max"],
        antibiotic_framework["psi_min"],
        antibiotic_framework["k"],
        antibiotic_framework["time frame"],
        antibiotic_framework["Initial Populations"]
        )

        fig, ax = plt.subplots(figsize=(14,6)) # prepare plot     

        # plot dynamics
        for initial_population in initial_populations:
            for name, params in genotypes.items():

                X = odeint(dX_dt, initial_population, time_frame,
                           args=(psi_max, psi_min, zMIC, k, params, self, antibiotic_concentration)) # args will be passed down to dX_dt
                ax.plot(time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')

        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e9)   

        # this is functionality for getting the legend out of the plot 
        # i dont use it because if i change the plot width now and 
        # try to plot env.variation later the plots will have different size
        # one option is to conditionally do it only if the function is called directly!

        # pos = ax.get_position() #returns bbox in order to manipulate width/height
        # ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        # ax.legend(bbox_to_anchor=(1.41, 1), fontsize="7") # place legend out of plot

        # only save when called directly
        if not is_called_from_another_function():
            save(f'./report/Population Dynamics', close=False)

        return fig, ax



        # Unpacking the dictionary into variables
        zMIC, antibiotic_concentration, psi_max, psi_min, k, time_frame, initial_populations = (
        antibiotic_framework["zMIC"],
        antibiotic_framework["Antibiotic Concentration"],
        antibiotic_framework["psi_max"],
        antibiotic_framework["psi_min"],
        antibiotic_framework["k"],
        antibiotic_framework["time frame"],
        antibiotic_framework["Initial Populations"]
        )

        fig, ax = plt.subplots(figsize=(14,6)) # prepare plot     

        # plot dynamics
        for initial_population in initial_populations:
            for name, params in genotypes.items():

                X = odeint(dENV_dt, [initial_population,0,0], time_frame,
                           args=(psi_max, psi_min, zMIC, k, params, self, antibiotic_concentration)) # args will be passed down to dX_dt
                ax.plot(time_frame, X[:,0], label=f'X0={'{:.0e}'.format(initial_population)} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')

        ax.set_xlabel('Time')
        ax.set_ylabel('Bacterial Density')
        ax.set_yscale('log')
        ax.set_ylim(1, 1e10)   

        antibiotic_exposure_layers_applier(time_frame,ax)

        # this is functionality for getting the legend out of the plot 
        # i dont use it because if i change the plot width now and 
        # try to plot env.variation later the plots will have different size
        # one option is to conditionally do it only if the function is called directly!

        # pos = ax.get_position() #returns bbox in order to manipulate width/height
        # ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
        # ax.legend(bbox_to_anchor=(1.41, 1), fontsize="7") # place legend out of plot

        # only save when called directly
        if not is_called_from_another_function():
            save(f'./report/New Population Dynamics', close=False)

        env = [dENV_dt(y, time, psi_max, psi_min, zMIC, k, params, self, antibiotic_concentration)[1] for time, y in zip(time_frame, X)]
        environmental_variation_layer_applier(time_frame, ax, env)
        save(f'./report/Dynamics With True Variation')

        fig, ax = plt.subplots(figsize=(14,6)) # prepare plot   



        ax.plot(time_frame, env, linestyle="dashdot", color="purple", label="True Environmental variation")
        ax.legend()
        ax.grid()
        save(f'./report/True Environemtal variation')


        return fig, ax


    def population_dynamics_antibiotic_frames(self, genotypes, antibiotic_framework):
        
        fig, ax = self.population_dynamics(genotypes, antibiotic_framework)

        time_frame = antibiotic_framework["time frame"]

        # add antibiotic exposure information
        antibiotic_exposure_layers_applier(time_frame,ax)

        # only save when called directly
        if not is_called_from_another_function():
            save(f'./report/Population Dynamics.Antibiotic Layers',close=False)

        return fig, ax

    def population_dynamics_antibiotic_frames_env_variation(self, genotypes, antibiotic_framework):
        
        fig, ax = self.population_dynamics_antibiotic_frames(genotypes,antibiotic_framework)

        time_frame = antibiotic_framework["time frame"]

        # add environmental variation information
        environmental_variation_layer_applier(time_frame, ax, self.variation)
        save(f'./report/Population Dynamics.Antibiotic Layers.Variation',close= False)

        return fig, ax

    def run_simulation(self, genotypes, antibiotic_framework):
        
        self.save()
        self.gene_reaction_norms(genotypes)
        self.gene_responses(genotypes)
        self.population_dynamics_antibiotic_frames_env_variation(genotypes, antibiotic_framework)


class Simulator():
    '''
    Simulator class is designed to process environments & genotypes.
    It produces plots for environmental variation, genotype responses & norms and population dynamics
    based on some predefined equations.
    '''
    def __init__(self, environment_params, genotype_params, antibiotic_framework):

        self.environment_params = environment_params
        self.genotypes = genotype_params
        self.antibiotic_framework = antibiotic_framework

        self.environments = self._yield_environments() # immediatly create environments
        print(f"simulator contains {len(self.environments)} environments & {len(self.genotypes)} genotypes.")

        #  # code will not work with one environment, so i put this condition
        # if len(self.environments) == 0 or len(self.environments) == 1:
        #     print("please provide more than 1 environment")
        #     exit()
    
    def run(self):

        print("Simulation begins...")

        self.envs_plot = self.yield_environment_plots_with_antibiotic_frames()
        self.norms_plot = self.yield_reaction_norms()
        self.responses_plot = self.yield_phenotypic_responses()
        self.dynamics_plot = self.yield_population_dynamics_with_antibiotic_frames_env_variation()

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
        
        row_vectors, rows, columns = self._plot_layer_constructor()

        # this technique is based on matplots basic tutorial
        # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.suptitle(u"Environmental Variations\nE\u209C = A·sin(2πt/LR) + B·ε", fontsize = 30, y= 0.95)
        fig.text(0.5, 0.07, "Time (t)", va="center",  fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Environmental variation (E)", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        # case of 1 environment
        if len(self.environments) == 1:
            fig, axs  = self.environments[0].fig, self.environments[0].ax
            plt.figure(fig) # activate the current figure in order to save correctly
        
        # check if grid has two dimensions (+unique values for the another parameter )
        elif axs.ndim > 1: 
            for row, vector in enumerate(row_vectors): # iterate on the stack of rows
                for column, env in enumerate(vector): # iterate on the row itself
                    custom_plot(axs[row,column],env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}', ylim=(0,1))
                    axs[row,0].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                    axs[-1,column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")
        
        # check if only the parameter of interest has different values
        elif len(row_vectors)>1:    
            # set label value for the one dimension 
            axs[-1].set_xlabel(f"R ={row_vectors[0][0].R}", rotation="horizontal", fontsize=14, weight="bold")
            for row, vector in enumerate(row_vectors):
                env = vector[0]
                custom_plot(axs[row],env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}', ylim=(0,1))
                axs[row].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")

        # else the parameters of interest has only one value and some other parameters has more values
        else:
            # set label value for the one dimension 
            axs[0].set_ylabel(f"A ={row_vectors[0][0].A}", rotation="horizontal", fontsize=14, weight="bold")
            for column, env in enumerate(row_vectors[0]): 
                custom_plot(axs[column],env.t, env.variation, label=f'A={env.A}\nB={env.B}\nL={env.L}\nR={env.R}\ntrimmed={env.trimmed}', ylim=(0,1))
                axs[column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")

        save('./report/Stacked Environments', close=False)
        return fig, axs

    def yield_reaction_norms(self):

        row_vectors, rows, columns = self._plot_layer_constructor()       
        
        fig = plt.figure(figsize=(columns*8, rows*6)) # empty figure for template
        gs = fig.add_gridspec(rows, columns, hspace=0.1, wspace=0) # grid with dimensions & space between plots
        axs = gs.subplots(sharey=True) # sharing the same y range (i think based on the bigger value)
        fig.text(0.35, 0.07, "Environmental variation (E)",   fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Response (I)", rotation="vertical", va="center", fontsize=20) # put only 1 y label
        fig.suptitle(u"Reaction Norms\n I=I\u2080 + b·C", fontsize = 30, y= 0.95)


        # case of 1 environment
        if len(self.environments) == 1:
            fig = self.environments[0].gene_reaction_norms(self.genotypes)
            plt.figure(fig) # activate the current figure in order to save correctly

        # check if grid has two dimensions (+unique values for the another parameter )
        elif axs.ndim > 1: 
            for row, vector in enumerate(row_vectors):
                for column, env in enumerate(vector):
                        for name, params in self.genotypes.items():
                            I = reaction_norm(params["I0"], params["b"], env.variation)
                            custom_plot(axs[row,column], env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
       
       # check if only the parameter of interest has different values          
        elif len(row_vectors)>1: 
             for row, vector in enumerate(row_vectors):
                env = vector[0]
                for name, params in self.genotypes.items():
                    I = reaction_norm(params["I0"], params["b"], env.variation)
                    custom_plot(axs[row], env.variation, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title = f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                    
        
        else:
            for column, env in enumerate(row_vectors[0]): 
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


        # case of 1 environment
        if len(self.environments) == 1:
            fig = self.environments[0].gene_responses(self.genotypes)
            plt.figure(fig) # activate the current figure in order to save correctly

        # check if grid has two dimensions (+unique values for the another parameter )
        elif axs.ndim > 1: 
            for row, vector in enumerate(row_vectors):
                for column, env in enumerate(vector):
                    custom_plot(axs[row,column], env.t, env.variation,  label='Environmental Variation', linestyle="dashdot", color="purple")
                    for name, params in self.genotypes.items():
                        I = reaction_norm(params["I0"], params["b"], env.variation)
                        custom_plot(axs[row,column], env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title=f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")
        
        
       # check if only the parameter of interest has different values          
        elif len(row_vectors)>1: 
             for row, vector in enumerate(row_vectors):
                env = vector[0]
                custom_plot(axs[row], env.t, env.variation,  label='Environmental Variation', linestyle="dashdot", color="purple")
                for name, params in self.genotypes.items():
                    I = reaction_norm(params["I0"], params["b"], env.variation)
                    custom_plot(axs[row], env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title=f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                           
            
        else:
            for column, env in enumerate(row_vectors[0]): 
                custom_plot(axs[column], env.t, env.variation,  label='Environmental Variation', linestyle="dashdot", color="purple")
                for name, params in self.genotypes.items():
                    I = reaction_norm(params["I0"], params["b"], env.variation)
                    custom_plot(axs[column], env.t, I, label=f"{name}, IO={params["I0"]}, b={params["b"]}", legend_title=f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}")                                                      

        save('./report/Stacked Phenotypics Responses')
        return fig
    
    def yield_population_dynamics(self):

        antibiotic_framework = self.antibiotic_framework

        # Unpacking the dictionary into variables
        zMIC, antibiotic_concentration, psi_max, psi_min, k, time_frame, initial_populations = (
            antibiotic_framework["zMIC"],
            antibiotic_framework["Antibiotic Concentration"],
            antibiotic_framework["psi_max"],
            antibiotic_framework["psi_min"],
            antibiotic_framework["k"],
            antibiotic_framework["time frame"],
            antibiotic_framework["Initial Populations"]
        )

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
                        f"c={antibiotic_concentration} \n"
                        r"$\bf{" + "Growth  Rate  Modifier:" + "}$" + f"{get_function_body(growth_rate_modifier)} \n"
                        r"$\bf{" + "Death  Rate  Modifier:" + "}$" + f"{get_function_body(death_rate_modifier)}"
                        ),
                        fontsize=20)
        fig.text(0.5, 0.07, "Time (t)", fontsize=20) # put only 1 x label
        fig.text(0.05, 0.5, "Bacterial density", rotation="vertical", va="center", fontsize=20) # put only 1 y label

        # check for only 1 environment
        if len(self.environments) == 1:
            fig, axs = self.environments[0].population_dynamics(self.genotypes, self.antibiotic_framework)
            plt.figure(fig) # activate the current figure in order to save correctly       

        # check if grid has two dimensions (+unique values for the another parameter )
        elif axs.ndim > 1: 
            for row, vector in enumerate(row_vectors):
                for column, env in enumerate(vector):
                    
                    for initial_population in initial_populations:
                        for name, params in self.genotypes.items():
                            X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibiotic_concentration)) # args will be passed down to dX_dt
                            custom_plot(axs[row,column], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))
                            axs[row,0].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                            axs[-1,column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")                       
       
       # check if only the parameter of interest has different values          
        elif len(row_vectors)>1: 
             for row, vector in enumerate(row_vectors):
                env = vector[0]
                for initial_population in initial_populations:
                    for name, params in self.genotypes.items():
                        X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibiotic_concentration)) # args will be passed down to dX_dt
                        custom_plot(axs[row], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))
                        axs[row].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                        axs[-1].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")                

        # or another parameter has different values            
        else:
            for column, env in enumerate(row_vectors[0]): 
                         
                for initial_population in initial_populations:
                    for name, params in self.genotypes.items():
                        X = odeint(dX_dt, initial_population, time_frame, args=(psi_max, psi_min, zMIC, k, params, env, antibiotic_concentration)) # args will be passed down to dX_dt
                        custom_plot(axs[column], time_frame, X, label=f'X0={'{:.0e}'.format(initial_population)} Genotype Params: I0={params["I0"]}, b={params["b"]}', legend_title= f" Environment Parameters: A={env.A}, B={env.B}, L={env.L}, R={env.R}", ylim=(1,1e10), yscale=('log'))                                
                        axs[0].set_ylabel(f"A ={env.A}", rotation="horizontal", fontsize=14, weight="bold")
                        axs[column].set_xlabel(f"R ={env.R}", rotation="horizontal", fontsize=14, weight="bold")
        
        save('./report/Stacked Population Dynamics', close=False)
        return fig, axs
    
    def yield_environment_plots_with_antibiotic_frames(self):
        # this function works right but has problem with the custom save function & that is why i dont use it.
        # problem is that custom save() doesnt take this fig as current and it saves an irrelevant plot of the past
        # when tried to activate fig with plt.figure(fig) saving it causes a traceback error with save() (somethings broken with plt.close())
        
        # Unpacking the dictionary into variables
        time_frame = (self.antibiotic_framework["time frame"])

        fig, axs = self.yield_environment_plots()

        # check for only 1 environment
        if len(self.environments) == 1:
            antibiotic_exposure_layers_applier(time_frame,axs)
        
        # check for two dimensions
        elif axs.ndim > 1:
            for vector in axs:
                for ax in vector:
                    antibiotic_exposure_layers_applier(time_frame,ax)

        # one dimension only            
        else:
            for ax in axs:
                antibiotic_exposure_layers_applier(time_frame,ax)

        save('./report/Stacked Environments with antibiotic layers', close=False)
        return fig
        
    def yield_population_dynamics_with_antibiotic_frames(self):
        # this function works right but has problem with the custom save function & that is why i dont use it.
        # problem is that custom save() doesnt take this fig as current and it saves an irrelevant plot of the past
        # when tried to activate fig with plt.figure(fig) saving it causes a traceback error with save() (somethings broken with plt.close())

        time_frame= (self.antibiotic_framework["time frame"])
        fig, axs = self.yield_population_dynamics()

        # check for only 1 environment
        if len(self.environments) == 1:
            print("indeed")
            antibiotic_exposure_layers_applier(time_frame,axs)
            

        elif axs.ndim > 1:
            for vector in axs:
                for ax in vector:
                    antibiotic_exposure_layers_applier(time_frame,ax)
        else:
            for ax in axs:
                antibiotic_exposure_layers_applier(time_frame,ax)

        save("./report/Stacked Population Dynamics with Antibiotics Layers", close=False)
        return fig, axs

    def yield_population_dynamics_with_antibiotic_frames_env_variation(self):

        antibiotic_framework = self.antibiotic_framework

        # Unpacking the dictionary into variables
        zMIC, antibiotic_concentration, psi_max, psi_min, k, time_frame, initial_populations = (
            antibiotic_framework["zMIC"],
            antibiotic_framework["Antibiotic Concentration"],
            antibiotic_framework["psi_max"],
            antibiotic_framework["psi_min"],
            antibiotic_framework["k"],
            antibiotic_framework["time frame"],
            antibiotic_framework["Initial Populations"]
        )
        
        row_vectors, rows, columns = self._plot_layer_constructor()

        fig, axs = self.yield_population_dynamics_with_antibiotic_frames()

        # check for only 1 environment
        if len(self.environments) == 1:
            fig, axs = self.environments[0].population_dynamics_antibiotic_frames_env_variation(self.genotypes, self.antibiotic_framework)
            plt.figure(fig) # activate the current figure in order to save correctly                           

        # check if grid has two dimensions (+unique values for another parameter )
        elif axs.ndim > 1:
            for row, vector in enumerate(row_vectors):
                for column, env in enumerate(vector):
                    
                    # add environmental variation information
                    environmental_variation_layer_applier(time_frame,axs[row,column],env.variation)

        # check if the the only parameter of interest has more than one value
        elif len(row_vectors)>1:
            for row, vector in enumerate(row_vectors):
                env = vector[0]

                # add environmental variation information
                environmental_variation_layer_applier(time_frame,axs[row],env.variation)
                         
        # else another parameter has variation            
        else:    
            for column, env in enumerate(row_vectors[0]): 

                # add environmental variation information
                environmental_variation_layer_applier(time_frame,axs[column],env.variation)
                                

        save('./report/Stacked Population Dynamics & antibiotic frames & env variation')
        return fig, axs


