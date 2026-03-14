import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from bokeh.plotting import figure, show, curdoc
from bokeh.models import Div, Spinner, TextInput
from bokeh.events import ValueSubmit, ButtonClick
from bokeh.layouts import layout, column, row

curdoc().theme = "dark_minimal"

def environmental_variation(params,t, epsilon):
    A, B, L, R = params


    return A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon

A = 1
B = 0
L = 1 # lifespan of ants
R = 1 # relative time of variation
T = 100

def plot(A=A, B=B, L=L, R=R, T=T):
        
        params = A, B, L, R 
        time = np.arange(T)
        time = np.linspace(0, T, T+1000)
        epsilon = np.random.normal(0, 1, T+1000)
        E = environmental_variation(params, time,epsilon)  

        p = figure(title = "Environmental variation",
                x_axis_label = "time", 
                y_axis_label = "Environment input value")

        p.width = 1280
        p.height = 720
        line = p.line(time, E, legend_label="Environmental variation function", line_width=3, color="red")

        return p

# widgets
time_spinner = Spinner(title="Time Duration", low=10, high=1000, step=1, value=T, width=200)
A_spinner = Spinner(title="A magnitude", low=0, high=1, step=0.1, value=A, width=200)
B_spinner = Spinner(title="B magnitude", low=0, high=1, step=0.1, value=B, width=200)
L_spinner = Spinner(title="L Lifespan (step per generation)", low=1, high=100, step=1, value=L, width=200)
R_spinner = Spinner(title="R Relevant variation frequency", low=1, high=100, step=1, value=R, width=200)

# Update functions
def update_plot(att, old, new):
    dashboard.children[1] = plot(A_spinner.value, B_spinner.value, L_spinner.value, R_spinner.value, time_spinner.value)

    return None

A_spinner.on_change("value", update_plot)
B_spinner.on_change("value", update_plot)
L_spinner.on_change("value", update_plot)
R_spinner.on_change("value", update_plot)
time_spinner.on_change("value", update_plot)
      


# layout = layout([[time_spinner, A_spinner, B_spinner, L_spinner,R_spinner],[plot()]])
widgets = row(time_spinner, A_spinner, B_spinner, L_spinner,R_spinner)
dashboard = column(widgets,plot())
curdoc().add_root(dashboard)









    
