import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import io 
import os
from PIL import Image 
import inspect
import ast
from matplotlib.lines import Line2D


from .equations import is_time_for_administration
# ╔══════════════════════════════════════════════════╗
# ║                    Tools                         ║
# ╚══════════════════════════════════════════════════╝

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
        ax.axvspan(i, i+1, facecolor='lightcoral', edgecolor='none', alpha=0.15 ) 
    for i in set(antibiotic_NOT_exposure_frame):
        ax.axvspan(i, i+1, facecolor='palegreen', edgecolor='none', alpha=0.15 )

    #create color patches for the legend to show
    exposure_patch = mpatches.Patch(color='red',  alpha=.2, label='Antibiotic Exposure')
    no_exposure_patch = mpatches.Patch(color='green', alpha=.2, label='No exposure')

    # https://www.statology.org/matplotlib-manual-legend/
    handles, labels = ax.get_legend_handles_labels() # extracting the previous legend stuff
    handles.extend([exposure_patch,no_exposure_patch]) # adding the patch to the old stuff
    ax.legend(handles=handles, bbox_to_anchor=(1.34, 1))

    return ax

def environmental_variation_layer_applier(time_frame, ax, variation):

    variation_axe = ax.twinx()

    # PLACE LEGEND OUT OF PLOT
    pos = variation_axe.get_position() #returns bbox in order to manipulate width/height
    variation_axe.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height]) # shrink figure's width in order to place legend outside of plot
    # variation_axe.legend(bbox_to_anchor=(1.34, 1)) # place legend out of plot
    
    # https://www.statology.org/matplotlib-manual-legend/
    #create color patches for the legend to show
    exposure_patch = mpatches.Patch(color='red',  alpha=.2, label='Antibiotic Exposure')
    no_exposure_patch = mpatches.Patch(color='green', alpha=.2, label='No exposure')
    variation_line = Line2D([0], [0], label='Nutrients Availability', color='black', alpha=0.2, linestyle='dashdot')
    handles, labels = ax.get_legend_handles_labels() # extracting the previous legend stuff
    handles.extend([exposure_patch, no_exposure_patch,variation_line]) # adding the patch to the old stuff
    ax.legend(handles=handles, bbox_to_anchor=(1.36, 1))
    
    custom_plot(variation_axe, time_frame, variation, linestyle="dashdot", color="black", alpha=0.2, ylim=(0,1))
    # variation_axe.yaxis.set_major_locator(ticker.NullLocator()) # remove ticks and labels rom y axis

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

def is_called_from_another_function():

    stack = inspect.stack()
    if len(stack) > 2:
        calling_function = stack[2].function
        if calling_function != "<module>":  # Exclude calls from module level
            return True
    return False

def generate_color_list(num_colors):
    # Get the default color cycle from Matplotlib
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    # If num_colors is greater than the length of the color cycle, repeat the cycle
    color_list = [color_cycle[i % len(color_cycle)] for i in range(num_colors)]
    
    return color_list