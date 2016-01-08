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
import struct
import argparse
import pickle
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

__author__ = "Christopher Hann-Soden"
__copyright__ = "Copyright 2015, Christopher Hann-Soden"
__credits__ = ["Christopher Hann-Soden"]
__licence__ = "GPL"
__version__ = "0.91"
__maintainer__ = "Christopher Hann-Soden"
__email__ = "channsoden@berkeley.edu"
__status__ = "Development"


def parse_arguments():
    global args
    parser = argparse.ArgumentParser(description='Generates ART read profiles from Illumina InterOp data sets. '\
                                     'Also creates linegraphs of quality profiles.')
    parser.add_argument('-p', '--profiles', action='store_true', help='only generate ART read profiles')
    parser.add_argument('-g', '--graphs', action='store_true', help='only generate line graphs of quality profiles')
    parser.add_argument('interopDirs', metavar='dirs', type=str, nargs='+', help='the list of input directories')
    parser.add_argument('-N', '--percentile', type=int, default=10, help='Nth percentile of quality scores to plot in the line graph')
    args = parser.parse_args()

    if args.profiles and args.graphs:
        sys.exit('Options -p and -g are incompatible with one another. Both profiles and graphs are produced by default.')

def main(interopDirs):
    check_for_required_files(interopDirs)
    interOp_data = {path: interop_quals(path) for path in interopDirs}

    if not args.profiles:
        # Calculate the mean and 10th percentile Q-scores at each cycle, then plot them in a line graph.
        lineseries_R1 = {}
        lineseries_R2 = {}
        for path, interop in interOp_data.items():
            lineseries_R1[path] = interop_meanNth_Q_by_position(interop[0], args.percentile)
            if len(interop) == 2:
                lineseries_R2[path] = interop_meanNth_Q_by_position(interop[1], args.percentile)
        linegraph(lineseries_R1, 'InterOp_quality_profile_R1.pdf')
        if lineseries_R2:
            linegraph(lineseries_R2, 'InterOp_quality_profile_R2.pdf')

    if not args.graphs:
        # Generate ART read profiles for each InterOp dataset.
        for path, interop in interOp_data.items():
            if len(interop) == 2:
                art_profile(interop[0], path+path[:-1]+'_R1.profile')
                art_profile(interop[1], path+path[:-1]+'_R2.profile')
            else:
                art_profile(interop[0], path+path[:-1]+'.profile')

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

def get_cycleinfo(RI_file):
    # Returns the total number of cycles in the run, the readlength, and a list containing the starting cycle numbers for the reads.
    tree = ET.parse(RI_file)
    root = tree.getroot()
    reads = root.iter('Read')
    
    cycles = 0
    read_idxs = []
    for read in reads:
        NumCycles = int(read.attrib['NumCycles'])
        if read.attrib['IsIndexedRead'] == 'N':
            readlength = NumCycles
            read_idxs.append(cycles)
        cycles += NumCycles

    return cycles, readlength, read_idxs

def parse_QMetrics(QM_file, cycles):
    # Returns distributions of quality scores by cycle number
    qmfh = open(QM_file)

    file_version, record_length, binning = struct.unpack('BBB', qmfh.read(3))

    if not file_version in [5, 6]:
        sys.exit('This QMetricsOut.bin is version {}.\nOnly versions 5 and 6 are supported.'.format(file_version))
    
    if binning:
        bins = struct.unpack('B', qmfh.read(1))[0]
        lower_bounds = list(struct.unpack('B'*bins, qmfh.read(bins)))
        upper_bounds = list(struct.unpack('B'*bins, qmfh.read(bins)))
        remapped_scores = list(struct.unpack('B'*bins, qmfh.read(bins)))
    else:
        bins = 50
        lower_bounds = range(1,50)
        upper_bounds = range(1,50)
        remapped_scores = range(1,51)

    Qdist = [np.zeros(bins) for x in range(cycles)]

    record = qmfh.read(record_length)
    while record:
        lane, tile, cycle = struct.unpack('HHH', record[:6])

        if file_version == 5:
            Qscore_counts = struct.unpack('I'*50, record[6:])
            binned_counts = np.array([Qscore_counts[score-1] for score in remapped_scores])
        elif file_version == 6:
            binned_counts = np.array(struct.unpack('i'*bins, record[6:]))

        Qdist[cycle-1] += binned_counts

        record = qmfh.read(record_length)

    for i, dist in enumerate(Qdist):
        if list(dist) == [0] * bins:
            # There exists a cycle for which no quality scores are recorded. This isn't right.
            print 'cycle', i+1
            sys.exit('Invalid Data: The file {} is missing data for cycle {} - aborting.'.format(QM_file, i+1))

    # Set the remapped scores as keys for the bins
    # Qdist[4][2] is the second bin in cycle 4
    # if remapped_scores == [7, 12, 17, 22, 27, 32, 37, 41]
    # then the second bin represents the counts of Q-score == 12
    Qdist = [dict(zip(remapped_scores, dist)) for dist in Qdist]

    return Qdist

def interop_quals(path):
    # Parses InterOp binary dataset, returns an ordered list of histogram-like distributions of quality scores for each cycle of the run.
    cycles, readlength, read_starts = get_cycleinfo(path+'RunInfo.xml')
 
    # Make an ordered list containing the frequency distribution of quality scores for each cycle
    qualdists = parse_QMetrics(path+'InterOp/QMetricsOut.bin', cycles)

    # Seperate out the reads
    read_qualdists = [qualdists[start:start+readlength] for start in read_starts]

    return read_qualdists

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
