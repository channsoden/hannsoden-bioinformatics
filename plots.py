#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

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


        
