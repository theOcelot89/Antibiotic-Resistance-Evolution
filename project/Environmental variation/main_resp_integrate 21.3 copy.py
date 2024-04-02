# ╔══════════════════════════════════════════════════╗
# ║                  Imports                         ║
# ╚══════════════════════════════════════════════════╝
#region
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Classes                         ║
# ╚══════════════════════════════════════════════════╝
#region

class Environment():

    def __init__(self, A = 1, B = 0.1, t = 365, L = 10, R = 10):
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

        self.variation[self.variation < 0] = 0  # if any of the values are negative, set them to 0
        self.variation[self.variation > 1] = 1    # if any of the values are greater than 1, set them to 1
        self.trimmed = True

        # reconstruct variation plot
        self.fig, self.ax = self._create_plot()

    def view(self):

        renderPeriod = 4
        initial_title = self.ax.get_title() # keep original title to put after rendering
        self.ax.set_title(self.ax.get_title() + "\n" + f'Will disappear after {renderPeriod} seconds') # temp title for rendering
        plt.show(block=False)
        plt.pause(renderPeriod)
        self.ax.set_title(initial_title) # set initial title again
        
    def save(self):
        self.fig.savefig(f"Env. variation t={len(self.t)}, A={self.A}, B={self.B}, L={self.L}, R={self.R}, trimmed={self.trimmed}.png")

#endregion

# ╔══════════════════════════════════════════════════╗
# ║                  Functions                       ║
# ╚══════════════════════════════════════════════════╝
#region
        
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

def dX_dt(X, t, psi_max, psi_min, zMIC, k, environment):
    '''function in which growth rate is calculated depending on the environmental conditions'''

    # decide in which timestep(e.g day) to quit the administration of antibiotic
    if t > 5: 
        a_t = 3 # antibiotic concentration
    else:
        a_t = 0
    
    current_env = environment.variation[int(t) % len(environment.t)] # E is the environmental variation (as an environmental Cue) at time t
    growth_rate_modifier = psi_max * reaction_norm(params["I0"], params["b"], current_env) # new psimax depending on plasticity
    growth_rate = np.log(10) * psi(a_t, growth_rate_modifier, psi_min, zMIC, k) * X

    return max(growth_rate, -X / 0.04)


#endregion

# ╔══════════════════════════════════════════════════╗
# ║                 Parameters                       ║
# ╚══════════════════════════════════════════════════╝
#region Genotypes

genotypes = {
    "Genotype 1": {"I0": 0, "b": 0},
    "Genotype 2": {"I0": 0.3, "b":0.5},
    "Genotype 3": {"I0": 0, "b": 1}
}
colors = ['blue', 'green', 'red'] # Define colors for each genotype

#endregion

#region antibiotic response curve

psi_min = -2 # maximum death rate
zMIC = 2 # concentration in which net growth rate is zero
k = 0.8  # Using a single mean k value
psi_max = 0.8  # maximal growth rate

# Time vector
t = np.linspace(0, 10, 10)
initial_populations = [1e3,1e4,1e5,1e6]

#endregion

# ╔══════════════════════════════════════════════════╗
# ║              Plots/Simulations                   ║
# ╚══════════════════════════════════════════════════╝

#region environment construction

env = Environment(A=1, B=0.1, L=10, R=100, t=110)
env.view()
env.save()
env.trim()
env.view()
env.save()

#endregion

#region responses to environmental variation

for (name, params), color in zip(genotypes.items(), colors):
    I = reaction_norm(params["I0"], params["b"], env.variation)
    env.ax.plot(env.t, I, label=f"{name},IO={params["I0"]},b ={params["b"]}", color=color)

env.ax.set_title('Phenotypic Response')
env.ax.set_xlabel('Time (t)')
env.ax.set_ylabel('Phenotypic response (I)')
env.ax.legend(bbox_to_anchor=(1.32, 1))

# Save the plot to a file
env.fig.savefig(f'Reaction Norms Across Time.png')

#endregion

#region reaction norms
    
fig, ax = plt.subplots(figsize=(12,6))

for (name, params), color in zip(genotypes.items(), colors):
    I = reaction_norm(params["I0"], params["b"], env.variation)
    plt.plot(env.variation, I, label=f"{name},IO={params["I0"]},b ={params["b"]}", color=color)

pos = ax.get_position() #returns bbox in order to manipulate width/height
ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 

ax.set_title('Reaction Norms for Genotypes')
ax.set_xlabel('Environmental Variation (E)')
ax.set_ylabel('Phenotypic Expression (I)')
# ax.legend()

fig.savefig('Reaction Norms for Genotypes.png') # Save the plot


#endregion

#region bacterial growth simulations

    #region different initial populations

fig, ax = plt.subplots(figsize=(14,6))

for X0 in initial_populations:

    X = odeint(dX_dt, X0, t, args=(psi_max, psi_min, zMIC, k, env)) # args will be passed down to dX_dt
    
    ax.plot(t, X, label=f'X0={'{:.0e}'.format(X0)} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')
    ax.set_xlabel('Time')
    ax.set_ylabel('Bacterial Density')
    ax.set_yscale('log')
    ax.set_ylim(1, 1e9)

pos = ax.get_position() #returns bbox in order to manipulate width/height
ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
ax.legend(bbox_to_anchor=(1.4, 1), fontsize="7") # place legend out of plot

fig.savefig(f' Different initial population dynamics.png')

    #endregion

    #region different genotypes

fig, ax = plt.subplots(figsize=(14,6))

for (name, params), color in zip(genotypes.items(), colors):
    
    X = odeint(dX_dt, initial_populations[0], t, args=(psi_max, psi_min, zMIC, k, env)) # args will be passed down to dX_dt
    
    ax.plot(t, X, label=f'X0={'{:.0e}'.format(initial_populations[0])} k={k}, Ψmax={psi_max}, Ψmin={psi_min}, MIC={zMIC}, I0={params["I0"]}, b={params["b"]} ')
    ax.set_xlabel('Time')
    ax.set_ylabel('Bacterial Density')
    ax.set_yscale('log')
    ax.set_ylim(1, 1e9)

pos = ax.get_position() #returns bbox in order to manipulate width/height
ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
ax.legend(bbox_to_anchor=(1.4, 1), fontsize="7") # place legend out of plot

fig.savefig(f' Different genotypes dynamics.png')
    #endregion

#endregion