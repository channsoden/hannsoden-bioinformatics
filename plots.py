#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection

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

def pretty_bar(ax, data, labels, title=None, shift = 0, barwidth=0.5,
               barcolor='gray', horizontal=False):
    if horizontal:
        plotter = ax.barh
        xlim = ax.set_ylim
        xticks = ax.set_yticks
        yticks = ax.set_xticks
        xticklabels = ax.set_yticklabels
        grid = ax.xaxis.grid
        rotation = -90
    else:
        plotter = ax.bar
        xlim = ax.set_xlim
        xticks = ax.set_xticks
        yticks = ax.set_yticks
        xticklabels = ax.set_xticklabels
        grid = ax.yaxis.grid
        rotation = 0

    # Matplotlib assumes you have numerical data for both axes
    # But for the X axis, we have categorical data
    # So we need to plug in a range of numbers to place the bars at
    x = list(range(len(data)))
    x = [i + shift for i in x]

    # Bars will appear at 0, 1, 2, etc
    # So set the x-limits to include these values
    xlim(-barwidth, len(data)-barwidth)

    # Here we turn the bars grey and remove their edge border.
    # Also make the bars a little more narrow.
    bars = plotter(x, data, barwidth, align='center', color=barcolor, edgecolor='none')

    # Set the ticks to match the bars and label them.
    if max([len(l) for l in labels]) > 4:
        rotation += 90
    xticks(x)
    xticklabels(labels, rotation=rotation)

    # Hide the frame around the plot
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    # Turn off the ticks
    ax.tick_params(bottom='off', top='off', left='off', right='off')
    # Overlay a white grid on the y axis
    grid(True, color='white', linestyle='solid')

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
             marker_alpha=0.3,
             plot_regression=True):

        with np.errstate(divide='ignore'):
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

        if plot_regression:
            if self.p_val <= self.alpha:
                line = '{}{}'.format(sig_style, sig_color)
            else:
                line = '{}{}'.format(insig_style, insig_color)
            
            regressx = np.linspace(self.x.min(), self.x.max(), num=100)
            prediction = regressx * self.slope + self.intercept
            with np.errstate(divide='ignore'):
                if logx:
                    regressx = np.log10(regressx)
                if logy:
                    prediction = np.log10(prediction)
            ax.plot(regressx, prediction, line)
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

def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
    # rotate theta such that the first axis is at the top
    theta += np.pi/2

    def draw_poly_patch(self):
        verts = unit_poly_verts(theta)
        return plt.Polygon(verts, closed=True, edgecolor='k')

    def draw_circle_patch(self):
        # unit circle centered on (0.5, 0.5)
        return plt.Circle((0.5, 0.5), 0.5)

    patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    if frame not in patch_dict:
        raise ValueError('unknown value for `frame`: %s' % frame)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = patch_dict[frame]

        def fill(self, *args, **kwargs):
            """Override fill so that line is closed by default"""
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            return self.draw_patch()

        def _gen_axes_spines(self):
            if frame == 'circle':
                return PolarAxes._gen_axes_spines(self)
            # The following is a hack to get the spines (i.e. the axes frame)
            # to draw correctly for a polygon frame.

            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)

            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}

    register_projection(RadarAxes)
    return theta


def unit_poly_verts(theta):
    """Return vertices of polygon for subplot axes.

    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts

def radar_plot(titles, data):
    theta = radar_factory(len(titles), frame='circle')
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='radar')

    ax.set_rgrids([0.2*i for i in range(1,6)])
    ax.set_ylim(0, 1)
    ax.set_title('this is a radar plot')
    ax.set_varlabels(titles)

    #ax.spines['polar'].set_linewidth(0) # to turn off bounding box
    [line.set_visible(False) for line in ax.xaxis.get_gridlines()]
    
    inspect = ax.xaxis.get_gridlines()
    print(inspect)
    print(dir(inspect))
    
    for d in data:
        ax.plot(theta, d)

    fig.savefig('radar_plot.png', bbox_inches='tight')
