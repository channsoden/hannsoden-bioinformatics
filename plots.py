#!/usr/bin/env python

from matplotlib import pyplot as plt

import numpy as np
import statsmodels.api as sm

import plotting_tools as pt

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

def pretty_bar(ax, data, labels, title=None, shift = 0, barwidth=0.5, barcolor='gray'):
    # Matplotlib assumes you have numerical data for both axes
    # But for the X axis, we have categorical data
    # So we need to plug in a range of numbers to place the bars at
    x = list(range(len(data)))
    x = [i + shift for i in x]

    # Bars will appear at 0, 1, 2, etc
    # So set the x-limits to include these values
    ax.set_xlim(-barwidth, len(data)-barwidth)

    # Here we turn the bars grey and remove their edge border.
    # Also make the bars a little more narrow.
    bars = ax.bar(x, data, align='center', color=barcolor, edgecolor='none', width=barwidth)

    # Set the ticks to match the bars and label them.
    if max([len(l) for l in labels]) > 4:
        rotation = 90
    else:
        rotation = 0
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=rotation)

    # Hide the frame around the plot
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    # Turn off the ticks
    ax.tick_params(bottom='off', top='off', left='off', right='off')
    # Overlay a white grid on the y axis
    ax.yaxis.grid(True, color='white', linestyle='solid')

    if title:
        ax.set_title(title, fontweight='bold')

    return bars

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
    yticks = list(range(max_y + 1))
    ax.set_yticks(yticks)

    xtick_labels = sorted(list(scaffold_positions.keys()), key=lambda k: scaffold_positions[k])
    xticks = [scaffold_positions[k] for k in xtick_labels]
    label_locations = [(xticks[i]+xticks[i+1])/2. for i in range(len(xticks)-1)]
    ax.set_xlim(0, max(xticks))
    ax.set_xticks(xticks)
    ax.set_xticklabels('') # turn off major tick labels
    ax.set_xticks(label_locations, minor=True)
    ax.set_xticklabels(xtick_labels, rotation='vertical', minor=True)
    ax.tick_params(axis='x', which='minor', bottom='off')

class regression_plot(object):
    def __init__(self, x, y, alpha = 0.05, label = None):
        self.x = x
        self.y = y
        self.label = label
        self.alpha = alpha

    def regress(self, slope = 'nonzero'):
        X = sm.add_constant(self.x)
        self.model = sm.OLS(self.y, X, missing='drop')
        results = self.model.fit()
        self.r2 = results.rsquared
        self.intercept, self.slope = results.params[:2]
        self.raw_int_p = results.pvalues[0]

        if slope == 'nonzero':
            self.raw_slope_p = results.pvalues[1]
        elif slope == 'negative' or slope == 'positive':
            # convert to single-tailed test for appropriate slope
            pval = results.pvalues[1] / 2.
            if ((slope == 'negative' and self.slope > 0) or
                (slope == 'positive' and self.slope < 0)):
                pval = 1. - pval
            self.raw_slope_p = pval
        else:
            raise ValueError("slope argument should be 'nonzero', 'negative', or 'positive'.")

        self.p_val = self.raw_slope_p
        return self.raw_slope_p

    def draw(self, ax,
             logx = False, logy = False,
             xlim = None, ylim = None,
             scientific = True,
             sig_color = 'k', insig_color = 'r',
             sig_style = '--', insig_style = '--',
             fit_report_location = None,
             fit_report_ha = 'left', fit_report_va = 'bottom',
             marker_alpha=0.3):

        if logx:
            X = np.log10(self.x)
        else:
            X = self.x
        if logy:
            Y = np.log10(self.y)
        else:
            Y = self.y

        ax.scatter(X, Y,  c='k', alpha = marker_alpha, edgecolors='none')

        if xlim:
            ax.set_xlim(*xlim)
        if ylim:
            ax.set_ylim(*ylim)
        if scientific:
            pt.scientific(ax)

        if self.p_val <= self.alpha:
            line = '{}{}'.format(sig_style, sig_color)
        else:
            line = '{}{}'.format(insig_style, insig_color)

        prediction = self.x * self.slope + self.intercept
        if logy:
            prediction = np.log10(prediction)
        ax.plot(X, prediction, line)
        if fit_report_location:
            fit_report = 'R2 = {:.4f}\np = {:.3E}'.format(self.r2, self.p_val)
            left, right = ax.get_xlim()
            bottom, top = ax.get_ylim()
            location = (left + (right-left) * fit_report_location[0],
                        bottom + (top-bottom) * fit_report_location[1])
            ax.text(location[0], location[1],
                    fit_report,
                    ha=fit_report_ha, va=fit_report_va)

        ax.set_title(self.label)
