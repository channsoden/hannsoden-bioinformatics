#!/usr/bin/env python

# From http://colorbrewer2.org/
brewer = ['#a6cee3',
          '#1f78b4',
          '#b2df8a',
          '#33a02c',
          '#fb9a99',
          '#e31a1c',
          '#fdbf6f',
          '#ff7f00',
          '#cab2d6',
          '#6a3d9a',
          '#ffff99',
          '#b15928']

def minimal(ax):
    #ax.set_axis_off() is better for everything
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Turn off ticks
    ax.tick_params(axis=u'both', which=u'both',length=0)

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

def line_text(ax, c1, c2, s, size=12, color='black', linewidth='2'):
    x1, y1 = c1
    x2, y2 = c2

    ax.text((x1+x2)/2., (y1+y2)/2., s, size=size, ha='center', va='bottom')
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth)
