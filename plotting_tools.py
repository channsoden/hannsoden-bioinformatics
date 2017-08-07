#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np

# From http://colorbrewer2.org/
brewer = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
          '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
          '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']

gioti = {'red': '#EF2636',
         'blue': '#3972C2',
         'green': '#2CE16C',
         'purple': '#944C7D'}

# 26 colors from "A Colour Alphabet and the Limits of Colour Coding"
# https://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
# Maximal contrast for a white background
alphabet = np.array(['#f0a3ff', '#0075dc', '#993f00', '#4c005c',
                     '#191919', '#005c31', '#2bce48', '#ffcc99',
                     '#808080', '#94ffb5', '#8f7c00', '#9dcc00',
                     '#c20088', '#003380', '#ffa405', '#ffa8bb',
                     '#426600', '#ff0010', '#5ef1f2', '#00998f',
                     '#e0ff66', '#740aff', '#990000', '#ffff80',
                     '#ffff00', '#ff5005'])

def minimal(ax, labels=False, ticks=False):
    #ax.set_axis_off() is better for everything
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Turn off ticks
    if not ticks:
        ax.tick_params(axis=u'both', which=u'both',length=0)

    if not labels:
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

def spine_size(ax, size):
    [i.set_linewidth(size) for i in ax.spines.itervalues()]
    ax.tick_params(axis='both', which='both', labelsize=size*5)
    
def scientific(ax, y=True, x=True):
    if y and not x:
        ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='y')
    elif x and not y:
        ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='x')
    elif x and y:
        ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='both')
    else:
        pass

def line_text(ax, c1, c2, s,
              size=12, color='black', linewidth='2', va='bottom',
              clip=True):
    x1, y1 = c1
    x2, y2 = c2

    ax.text((x1+x2)/2., (y1+y2)/2., s, size=size, ha='center', va=va)
    line = plt.Line2D([x1, x2], [y1, y2], color=color, linewidth=linewidth)
    line.set_clip_on(clip)
    ax.add_line(line)
    #ax.plot()

def color_legend(ax, labels, colors, position='upper right'):
    import matplotlib.patches as mpatches
    
    handles = [mpatches.Patch(color=c, label=l) for l, c in zip(labels, colors)]
    ax.legend(handles=handles, loc=position, frameon=False)

def color_scale(value, max_color=np.array([0., 0., 0.]), min_color=np.array([1., 1., 1.])):
    dif = max_color - min_color
    color = min_color + value * dif
    return color

def heat_scale(array):
    # 0 = blue = (0, 0, 1)
    # 0.5 = white = (1, 1, 1)
    # 1 = red = (1, 0, 0)
    colors = np.empty((len(array), 3))
    colors[:, 0] = 2 * array # Red
    #colors[:, 0][colors[:, 0] < 0] = 0
    colors[:, 1] = 1 - 2 * np.abs(array - 0.5) # Green
    colors[:, 2] = 2 * (1 - array) # Blue
    colors[colors > 1] = 1.
    return colors

def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def hex_to_rgb(hex_color):
    hex_color = hex_color.lower()
    hex_color = hex_color.strip('#').split('x')[-1]
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return r, g, b
