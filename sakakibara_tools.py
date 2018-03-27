#!/usr/bin/env python
# Standard modules
import os, bisect
import subprocess as sp

# Nonstandard modules
import pandas as pd

# My modules
from fasta_tools import get_absolute_positions

def murasaki_pairwise(pair, outdir, outprefix, weight = 28, length = 36):
    anchors = '{}/{}.anchors'.format(outdir, outprefix)
    if not os.path.isfile(anchors):
        genomes = ' '.join(list(pair.values()))
        command = "murasaki -p[{}:{}] -d {} -n {} -H2 {}".format(weight, length, outdir, outprefix, genomes)
        # if it's taking forever try -m [int] to skip hashes that appear [int] times
        sp.Popen(command, shell=True).wait()
    else:
        print('Murasaki output already present at {}: using existing file.'.format(anchors))
    return anchors

def abs_to_rel(anchors, outfile, verbose = True):
    # Converts Murasaki output format using genome-based positions to OSFinder input format using chromosome/scaffold/contig-based positions
    stub = anchors[:-len('.anchors')]

    if not os.path.isfile(outfile):
        genome_file = stub + '.seqs'
        genomes, intab = read_sakakibara(anchors, genome_file, 'murasaki')
        out_header, sep = sakakibara_header(genomes, 'osfinder')
        outtab = pd.DataFrame(index = intab.index, columns=out_header)

        if verbose:
            print('Munging {} genomes from Murasaki format to OSFinder format.'.format(len(genomes)))

        for genome in genomes:
            genome_positions = sorted(get_absolute_positions(genome).values())
            scaf = genome+'.scaf'
            start = genome+'.start'
            stop = genome+'.stop'
            sign = genome+'.sign'
            
            # Convert Murasaki signed coordinates to OSFinder unsigned coordinates.
            # Murasaki position -1 is complementary to 1
            # -2 is complementary to 2, . . . -SIZE is comp to SIZE
            # -1 is NOT complementary to SIZE
            
            intab[start] = intab[start].apply(abs)
            intab[stop] = intab[stop].apply(abs)
            
            def order(row):
                # If the coord is negative need to swap start and stop
                # such that start < stop.
                if row[start] > row[stop]:
                    row[start], row[stop] = row[stop], row[start]
                return row
        
            intab = intab.apply(order, axis=1)
            
            #for index, row in intab.iterrows():
            #    if intab.loc[index, start] > intab.loc[index, stop]:
            #        intab.set_value(index, start, intab.loc[index, stop])
            #        intab.set_value(index, stop, intab.loc[index, start])

            # get number of scaffold and position adjustment for scaffold 
            outtab[scaf] = intab[start].apply(scaffold_number, args = (genome_positions, ))
            delta = outtab[scaf].apply(lambda x: genome_positions[x-1] - 1)
            
            # delta converts to relative position
            outtab[start] = intab[start] - delta
            outtab[stop] = intab[stop] - delta
            
            outtab[sign] = intab[sign]

            if verbose:
                print("\tFinished munging " + genome)
            
        outtab.to_csv(outfile, sep = sep, header = False, index = False)
        print('OSFinder format written to {}.'.format(outfile))
    else:
        print('OSFinder input found at {}: using existing file.'.format(outfile))

    return outfile

def read_sakakibara(file_name, genome_file, fmt):
    # Reads Murasaki and OSFinder formated coordinate files into pandas dataframes
    genomes = [line.strip() for line in open(genome_file, 'r')]
    header, sep = sakakibara_header(genomes, fmt)
    table = pd.read_csv(file_name, sep = sep, names = header, header = None, index_col = False)
    return genomes, table

def sakakibara_header(genomes, fmt):
    header = []
    if fmt == 'murasaki':
        [header.extend([genome+'.start', genome+'.stop', genome+'.sign']) for genome in genomes]
        sep = '\t'
    elif fmt == 'osfinder':
        [header.extend([genome+'.scaf', genome+'.start', genome+'.stop', genome+'.sign']) for genome in genomes]
        sep = ' '
    else:
        raise ValueError("fmt must be 'murasaki' or 'osfinder'")

    return header, sep

def scaffold_number(position, genome_positions):
    return bisect.bisect(genome_positions, position)

def osfinder_pairwise(anchors, outprefix, min_len = 1000):
    segments = outprefix + '.os.txt'

    if not os.path.isfile(segments):
        command = 'osfinder -i {} -o {} -n 2 -s {}'.format(anchors, outprefix, min_len)
        sp.Popen(command, shell=True).wait()
    else:
        print('OSFinder output already present at {}: using existing file.'.format(segments))
    return segments

def osfinder_to_grimm(osfinder_results, genome_file, outfile, genome_abbreviations):
    genomes, tab = read_sakakibara(osfinder_results, genome_file, 'osfinder')
    
    grimm_genomes = {genome: [] for genome in genomes}
    for genome in genomes:
        scaf = genome+'.scaf'
        start = genome+'.start'
        sign = genome+'.sign'

        tab = tab.sort_values([scaf, start])

        current_scaf = tab.iloc[0][scaf]
        for i, segment in tab.iterrows():
            if segment[scaf] != current_scaf:
                # End of the scaffold
                grimm_genomes[genome].append('$')
                current_scaf = segment[scaf]

            grimm_genomes[genome].append(str(eval(segment[sign]+str(i+1)))) # i+1 because GRIMM cannot start from 0

    outfh = open(outfile, 'w')
    for genome, segs in list(grimm_genomes.items()):
        outfh.write(' '.join(['>', genome_abbreviations[genome], '\n']))
        outfh.write(' '.join(segs) + '\n')
        outfh.write('\n')
    outfh.close()
