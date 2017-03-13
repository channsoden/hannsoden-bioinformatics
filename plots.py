#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import numpy as np

def box_plot(ax, data, color='black', spacing = 1, offset = 0, sym = 'k+'):
    categories = sorted(data.keys())
    series = [data[category] for category in categories]
    positions = [1+(i*spacing)+offset for i in range(len(series))]

    box = ax.boxplot(series, patch_artist=True, positions=positions, sym=sym, widths=0.4, vert=True)
    [patch.set_facecolor(color) for patch in box['boxes']]
    [patch.set_color(color) for patch in box['boxes']]
    [line.set_color('white') for line in box['medians']]
    [line.set_color(color) for line in box['whiskers']]
    [line.set_color(color) for line in box['caps']]
    [sym.set_color(color) for sym in box['fliers']]

    ax.set_xticklabels(categories)
        
def manhattan(ax, positions, p_vals, scaffold_positions, color = 'black', sig = 2):
    rgba = np.zeros((len(positions), 4))
    if type(color) == str:
        rgba[:, 0:3] = colors.to_rgba(color)[0:3]
    else:
        rgba[:, 0:3] = color
    #rgba[p_vals < sig, 0:3] = np.array([0., 0., 0.]) # insignificant points are black
    rgba[:, 3] = 1 # make opaque
    
    ax.scatter(positions, p_vals, color = rgba, s = 5)

    ax.set_ylabel('-log(p)')

    max_y = int(max(p_vals)) + 1
    ax.set_ylim(-0.3, max_y)
    yticks = range(max_y + 1)
    ax.set_yticks(yticks)

    xtick_labels = sorted(scaffold_positions.keys(), key=lambda k: scaffold_positions[k])
    xticks = [scaffold_positions[k] for k in xtick_labels]
    label_locations = [(xticks[i]+xticks[i+1])/2. for i in range(len(xticks)-1)]
    ax.set_xlim(0, max(xticks))
    ax.set_xticks(xticks)
    ax.set_xticklabels('') # turn off major tick labels
    ax.set_xticks(label_locations, minor=True)
    ax.set_xticklabels(xtick_labels, rotation='vertical', minor=True)
    ax.tick_params(axis='x', which='minor', bottom='off')


