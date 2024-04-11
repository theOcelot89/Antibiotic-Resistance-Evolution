import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

#define data to plot
x = [1, 2, 3, 4, 5, 6, 7]
y = [2, 3, 5, 8, 12, 18, 27]

#create scatter plot of x vs. y
plt.scatter(x, y, label='Original Data', color='steelblue')

#define handles and labels that will get added to legend
handles, labels = plt.gca().get_legend_handles_labels()

print(handles, labels)

#define patches and lines to add to legend
patch1 = mpatches.Patch(color='orange', label='First Manual Patch')
patch2 = mpatches.Patch(color='orange', label='First Manual Patch')   
line1 = Line2D([0], [0], label='First Manual Line', color='purple')
line2 = Line2D([0], [0], label='Second Manual Line', color='red')

#add handles
handles.extend([patch1, line1, line2])

#add legend
plt.legend(handles=handles)

#display plot
plt.show()