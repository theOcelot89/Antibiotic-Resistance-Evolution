import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint



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

        # yield variation
        self.variation = self.A * np.sin(2 * np.pi * self.t / (self.L * self.R)) + self.B * self.epsilon

        #construct plot 
        def plot(self):
            fig, ax = plt.subplots(figsize=(12, 6))   
            ax.plot(self.t, self.variation, label='Environmental Variation')
            pos = ax.get_position() #returns bbox in order to manipulate width/height
            ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
            ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 
            ax.set_title(f'A={self.A}, B={self.B}, L={self.L}, R={self.R} ')

            return fig

        self.plot = plot(self)

        print(f"New environment created!")

    def fixate(self):
        '''limiting environmental variation between 0 & 1'''

        self.variation[self.variation < 0] = 0  # if any of the values are negative, set them to 0
        self.variation[self.variation > 1] = 1    # if any of the values are greater than 1, set them to 1

    def view(self):
        plt.show(block=False)
        plt.pause(3)

    def save(self):
        self.plot.savefig(f"Env Variation A={self.A}, B={self.B}, L={self.L}, R={self.R}.png")

    # def view(self):
    #     fig, ax = plt.subplots(figsize=(12, 6))   
    #     ax.plot(self.t, self.variation, label='Environmental Variation')
    #     pos = ax.get_position() #returns bbox in order to manipulate width/height
    #     ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
    #     ax.legend(bbox_to_anchor=(1.32, 1)) # place legend out of plot 
    #     ax.set_title(f'A={self.A}, B={self.B}, L={self.L}, R={self.R} ')
    #     plt.show(block=False)
    #     plt.pause(3)

    #     return fig

        # fig.savefig(f'Environmental variation.png')



env = Environment()
env.view()
env.save()
x = None

# print(help(env))


