#!/usr/bin/env python
"""
A script to create ART read profiles and graph the Q-score profiles
from Illumina InterOp datasets.

./InterOp_profiler.py -h for usage information
    
The input datasets should each be in their own directories, and the
directory names should be passed to this script. Each directory must
contain a RunInfo.xml file and InterOp/ directory containing
QMetricsOut.bin for a different Illumina sequencing run.
"""
import sys
import os
import argparse
import pickle
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from illuminate import InteropDataset

__author__ = "Christopher Hann-Soden"
__copyright__ = "Copyright 2015, Christopher Hann-Soden"
__credits__ = ["Christopher Hann-Soden"]
__licence__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Christopher Hann-Soden"
__email__ = "channsoden@berkeley.edu"
__status__ = "Development"


def parse_arguments():
    global args
    parser = argparse.ArgumentParser(description='Generates ART read profiles from Illumina InterOp data sets. '\
                                     'Also creates linegraphs of quality profiles.')
    parser.add_argument('-p', '--profiles', action='store_true', help='only generate ART read profiles')
    parser.add_argument('-g', '--graphs', action='store_true', help='only generate line graphs of quality profiles')
    parser.add_argument('-f', '--forcerestart', action='store_true', help='forces reparsing of InterOp dataset')
    parser.add_argument('interopDirs', metavar='dir', type=str, nargs='+', help='the list of input directories')
    parser.add_argument('-N', '--percentile', type=int, default=10, help='Nth percentile of quality scores to plot in the line graph')
    args = parser.parse_args()

    if args.profiles and args.graphs:
        sys.exit('Options -p and -g are incompatible with one another. Both profiles and graphs are produced by default.')

def main(interopDirs):
    check_for_required_files(interopDirs)
    interOp_data = {path: parse_interop_quals(path) for path in interopDirs}

    if not args.profiles:
        # Calculate the mean and 10th percentile Q-scores at each cycle, then plot them in a line graph.
        lineseries_R1 = {}
        lineseries_R2 = {}
        for path, interop in interOp_data.items():
            qualdists, readlength, PE = interop
            if PE:
                qualdists_R1, qualdists_R2 = split_reads_interop_distribution(qualdists, readlength)
                lineseries_R1[path] = interop_meanNth_Q_by_position(qualdists_R1, args.percentile)
                lineseries_R2[path] = interop_meanNth_Q_by_position(qualdists_R2, args.percentile)
            else:
                qualdists = qualdists[:readlength]
                lineseries_R1[path] = interop_meanNth_Q_by_position(qualdists, args.percentile)
        #lineseries = {path: interop_meanNth_Q_by_position(interOp_data[path][0], 10) for path in interopDirs}
        linegraph(lineseries_R1, 'InterOp_quality_profile_R1.pdf')
        if lineseries_R2:
            linegraph(lineseries_R2, 'InterOp_quality_profile_R2.pdf')

    if not args.graphs:
        # Generate ART read profiles for each InterOp dataset.
        for path, interop in interOp_data.items():
            qualdists, readlength, PE = interop
            if PE:
                qualdists_R1, qualdists_R2 = split_reads_interop_distribution(qualdists, readlength)
                art_profile(qualdists_R1, path+path[:-1]+'_R1.profile')
                art_profile(qualdists_R2, path+path[:-1]+'_R2.profile')
            else:
                qualdists = qualdists[:readlength]
                art_profile(qualdists, path+path[:-1]+'.profile')

def data_present(path):
    # Should make sure InterOp/ and RunInfo.xml are present, and that QMetricsOut.bin is under InterOp/
    return (os.path.isdir(path+'InterOp/') and
            os.path.isfile(path+'RunInfo.xml') and
            os.path.isfile(path+'InterOp/QMetricsOut.bin')
            )

def check_for_required_files(directories):
    missing_dirs = [path for path in directories if not data_present(path)]
    if missing_dirs:
        sys.stderr.write('InterOp/QmetricsOut.bin or RunInfo.xml not found in the following input directories:\n')
        for missing_dir in missing_dirs:
            sys.stderr.write('\t'+missing_dir+'\n')
        sys.stderr.write('Please make sure all input directories contain these required files.')
        sys.exit(1)

def parse_interop_quals(path):
    # Parses InterOp binary dataset, returns an ordered list of histogram-like distributions of quality scores for each cycle of the run.
    if (os.path.isfile(path+'qualdists.pickle') and os.path.isfile(path+'readlength.pickle') and
        os.path.isfile(path+'PE.pickle') and not args.forcerestart):
        qualdists = pickle.load(open(path+'qualdists.pickle', 'r'))
        readlength = pickle.load(open(path+'readlength.pickle', 'r'))
        PE = pickle.load(open(path+'PE.pickle', 'r'))
    else:
        myDataset = InteropDataset(path)
        qualitymetrics = myDataset.QualityMetrics()

        readlength = qualitymetrics.read_tiers[0]
        PE = len(qualitymetrics.read_tiers) == 3
        
        cycles = set(qualitymetrics.df.cycle)
        # Make an ordered list containing the frequency distribution of quality scores for each cycle
        qualdists = []
        for cycle in cycles:
            qualdist = dict(qualitymetrics.df[(qualitymetrics.df.cycle == cycle)].sum()[:-3])
            qualdist = {int(key.strip('q')):value for key, value in qualdist.items()}
            qualdists.append(qualdist)
        pickle.dump(qualdists, open(path+'qualdists.pickle', 'w'))
        pickle.dump(readlength, open(path+'readlength.pickle', 'w'))
        pickle.dump(PE, open(path+'PE.pickle', 'w'))

    return (qualdists, readlength, PE)

def split_reads_interop_distribution(qualdists, readlength):
    # For paired end run only.
    # Returns two qualdists corresponding to the first and second reads.
    # The index is read in between R1 and R2, so this is the first N and last N cycles.
    return qualdists[:readlength], qualdists[-readlength:]

def interop_meanNth_Q_by_position(qualdists, N):
    # Return an ordered list of the mean and Nth percentile Q-scores for each cycle in an Illumina run from an InterOp dataset.
    means = [hist_mean(qualdist) for qualdist in qualdists]
    tenths = [hist_percentile(qualdist, N) for qualdist in qualdists]
    return zip(means, tenths)

def hist_mean(histogram):
    # Takes a histogram style dictionary, with some values as keys, and the frequencies of occurance of those values as values.
    # Returns the mean value of the distribution.
    stacks = [key * value for key, value in histogram.items()]
    N = float(sum(histogram.values()))
    total = float(sum(stacks))
    return (total / N)

def hist_percentile(histogram, percent):
    # Takes a histogram style dictionary, with some values as keys, and the frequencies of occurance of those values as values.
    # Returns Nth percentile of the distribution (i.e. N% of occurances are at less than this value).
    total = sum(histogram.values())
    Nth = total * percent / 100
    values = histogram.keys()
    values.sort()
    pos = 0
    for value in values:
        pos += histogram[value]
        if pos >= Nth:
            return value

def linegraph(series_dict, output_name):
    # Takes a dictionary of data series. Each data series consistes of a list of tuples. tuple[0] is the mean, tuple[1] is the Nth percentile.
    # series_dict = ['series1':[(X1, P1), (X2, P2), ...], 'series2':[ ... ], ...]
    # Different data series are plotted in different colors.
    # Mean values are plotted as solid lines. Percentile values are plotted as dotted lines.
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111)
    ax.set_ylim(0, 50)
    ax.set_xlim(0, len(max(series_dict.values())))
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Q-score')

    # Tableu Color Blind 10 scheme
    colors = [(255,128,14), (171,171,171), (95,158,209), (89,89,89), (0,107,164),
              (255,188,121), (207,207,207), (200,82,0), (162,200,236), (137,137,137)]
    colors = [map(lambda x: x/255., triplet) for triplet in colors]

    legend_styles = []
    legend_labels = []
    for i, label in enumerate(series_dict.keys()):
        i = i % 10
        color = colors[i]

        legend_styles.append(plt.Line2D((0,1),(0,0), color=color, linestyle='-'))
        legend_styles.append(plt.Line2D((0,1),(0,0), color=color, linestyle=':'))
        legend_labels.append(label[:-1]+' mean')
        legend_labels.append(label[:-1]+' {}th %'.format(str(args.percentile)))
        
        for x in range(len(series_dict[label][:-1])):
            ax.plot([x, x+1], [series_dict[label][x][0], series_dict[label][x+1][0]],
                    lw=1, color=color, linestyle='-', label=label+'_mean')
            ax.plot([x, x+1], [series_dict[label][x][1], series_dict[label][x+1][1]],
                    lw=1, color=color, linestyle=':', label=label+'_Nth_percentile')

    legend = ax.legend(legend_styles, legend_labels, loc='upper right')

    plt.savefig(output_name, bbox_inches='tight')

def art_profile(qualdist, outfile):
    # Creates ART sequencing profiles from the parsed InterOp quality distribution.
    # Because the InterOp datasets don't contain base calling information, I approximate the profile by assuming equal distribution
    # of the bases and no Ns. This is obviously not true, but unavoidable.
    bases = ['.', 'A', 'T', 'G', 'C']
    
    outfh = open(outfile, 'w')

    for base in bases:
        for i, cycle in enumerate(qualdist):
            reduced = {key: value for key, value in cycle.items() if value != 0}
            line1 = [base] + [str(i)] + map(str, sorted(reduced.keys()))
            line2 = [base] + [str(i)] + map(lambda x: str(x) if base == '.' else str(x/4), [cycle[Qscore] for Qscore in sorted(reduced.keys())])
            outfh.write('\t'.join(line1) + '\n')
            outfh.write('\t'.join(line2) + '\n')

    outfh.close()

if __name__ == '__main__':
    parse_arguments()
    main(args.interopDirs)
