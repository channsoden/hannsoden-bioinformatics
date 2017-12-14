#!/usr/bin/env python
import re, sys, time
from bisect import bisect_right

# hannsoden_bioinformatics
from processing_tools import mapPool

###############################################################
## Rewritten, fixed, and optimized by Christopher Hann-Soden ##
## https://github.com/channsoden/hannsoden-bioinformatics    ##
###############################################################

class char_set(object):
    def __init__(self, unknowns, alphabet):
        self.dna = {'A': 'ARWMDHVNarwmdhv',
                    'T': 'TUYWKBDHNtuywkbdhn',
                    'C': 'CYSMBHVNcysmbhvn',
                    'G': 'GRSKBDVNgrskbdvn',
                    'a': 'ARWMDHVNarwmdhv',
                    't': 'TUYWKBDHNtuywkbdhn',
                    'c': 'CYSMBHVNcysmbhvn',
                    'g': 'GRSKBDVNgrskbdvn',
                    '-':'-'}

        self.protein = {'A': 'Aa',  'C': 'Cc',  'D': 'Dd',  'E': 'Ee',
                        'F': 'Ff',  'G': 'Gg',  'H': 'Hh',  'I': 'Ii',
                        'K': 'Kk',  'L': 'Ll',  'M': 'Mm',  'N': 'Nn',
                        'P': 'Pp',  'Q': 'Qq',  'R': 'Rr',  'S': 'Ss',
                        'T': 'Tt',  'V': 'Vv',  'W': 'Ww',  'Y': 'Yy',
                        'a': 'Aa',  'c': 'Cc',  'd': 'Dd',  'e': 'Ee',
                        'f': 'Ff',  'g': 'Gg',  'h': 'Hh',  'i': 'Ii',
                        'k': 'Kk',  'l': 'Ll',  'm': 'Mm',  'n': 'Nn',
                        'p': 'Pp',  'q': 'Qq',  'r': 'Rr',  's': 'Ss',
                        't': 'Tt',  'v': 'Vv',  'w': 'Ww',  'y': 'Yy',
                        '-': '-'}

        self.standard = {chr(i):chr(i) for i in range(256)}

        self.dna = {k:v+unknowns for k,v in self.dna.items()}
        self.protein = {k:v+unknowns for k,v in self.protein.items()}
        self.standard = {k:v+unknowns for k,v in self.standard.items()}

        self.name = alphabet
        self.alphabet = eval('self.'+alphabet)

        if alphabet == 'dna':
            unknowns += 'RYSWKMBDHVNryswkmbdhvn'
        self.unknown = set(list(unknowns))

    def __repr__(self):
        return self.name

class bitsize(object):
    # This class's properties are global
    pass

def set_bitsize(n):
    bitsize.n = n
    bitsize.max = (1 << n) -1

def report(message):
    now = time.localtime()
    print '{}:{}:{} - {}'.format(now.tm_hour, now.tm_min, now.tm_sec, message)
    
def partition_sites(seqs, args):
    sites = [[s[i] for s in seqs] for i in range(len(seqs[0]))]
    jobs = [(partition, (site, args.alphabet)) for site in sites]
    partitions = mapPool(args.threads, jobs, daemonic=True, chunksize=10000)
    
    patterns = list( set(partitions) )
    partitions = [patterns.index(part) for part in partitions]
    return partitions, patterns

def partition(site, alphabet):
    # Unknown sites should not cause disagreement.
    # Unknown sites should thus be present in multiple
    # partitions, depending on what the unknown character is.
    traits = set(site) - alphabet.unknown
    part = [partcode(site, trait, alphabet) for trait in traits]
    # Catch sequences that aren't in any partition.
    # These are unknown or ambiguous characters that are
    # not the same as any known character in the column.
    leftovers = not_any(part)
    if leftovers:
        part.append(leftovers)
    part.sort()
    return tuple(part)

def partcode(site, trait, alphabet):
    # Returns a binary encoded subset of site where
    # 1 bits mean identiy and 0 bits mean nonidentity
    b = 0
    for c in site:
        b = b << 1
        if c in alphabet.alphabet[trait]:
            b +=1
    return b

def not_any(binaries):
    # Returns bits where none are 1 = 1 and all other bits = 0.
    # bitsize must be set, and all integers in binaries must be under bitsize.
    NE = 0
    for b in binaries:
        NE = NE | b

    return bitsize.max ^ NE

def calculate_rates(patterns, pattern_counts, nMinusOne, num_invariants, invariant_index, partitions, args):
    # parallelize, since this step can be very long
    jobs = [(score_conflict, (pat, patterns, pattern_counts, nMinusOne, num_invariants))
            for pat in patterns]
    pattern_conflicts = mapPool(args.threads, jobs, daemonic=True, chunksize=100)
    
    pattern_conflicts[invariant_index] = 0 # definitionally, and above calculation doesn't account for invariant sites
    pattern_rates = [1.-c for c in pattern_conflicts]
    # Expand pattern_rates into places where they occur in partitons.
    rates = [pattern_rates[i] for i in partitions]
    return pattern_rates, rates

def score_conflict(patA, patterns, pattern_counts, nMinusOne, num_invariants):
    # patterns are the unique partiton patterns
    # pattern_counts are the number of times those patterns appear in partitions
    # multiply each conflict score by the number of times it appears in partitions.
    score = -num_invariants
    for patB, count in zip(patterns, pattern_counts):
        conflictions = sum( [conflict(gb, patA) for gb in patB] )
        pat_score = conflictions / len(patB)
        score += pat_score * count

    return score / nMinusOne

def conflict(groupB, patA):
    # Returns 1. if groupB is not a subset of any group in patA
    # Returns 0. if groupB is a subset of a group in patA
    for groupA in patA:
        if binary_conflict(groupA, groupB):
            return 1.
    return 0.

def binary_conflict(x, y):
    # Returns true if x has some, but not all bits in y
    union = x & y
    return union and union != y

def linear_binning(rates, args):
    binsize = 1. / args.bins
    bins = [i * binsize for i in range(args.bins)]
    return hist(bins, rates)


def rate_partitioning(pattern_rates, pattern_counts, rates, args):
    upper = max(pattern_rates)
    lower = min(pattern_rates)
    d = args.bins
    assert d > 1, '"--bins" must be set greater than 1 when using "-bt rota"'

    rate_frequencies = sorted(zip(pattern_rates, pattern_counts), reverse=True)
    
    thresh = len(rates) * 0.9 # fastest bin contains <10% of all sites
    total = 0
    for rate, count in rate_frequencies:
        total += count
        if total > thresh:
            break
    # if I include rate in the lowest bin, it will be more than 10%
    # so I cannot include rate

    bins = [upper, upper - ((upper - lower) / d)]
    d += 0.3
    while bins[-1] >= rate:
        d += 0.3
        upper = upper - ((upper - lower) / d)
        bins.append(upper)
    bins.append(lower)
    bins = bins[-2::-1]

    return hist(bins, rates)    

def hist(bins, rates):
    hist = [0 for b in bins]
    ranks = []
    partitions = [[] for b in bins]
    for i, r in enumerate(rates):
        rank = bisect_right(bins, r) - 1
        hist[rank] += 1
        ranks.append(rank)
        partitions[rank].append(i)
    return bins, hist, ranks, partitions

def output_phylip(filename, names, seqs):
    # Uses relaxed phylip format
    fh = open(filename, 'w')

    # Print the header.
    n = len(seqs)
    m = len(seqs[0])
    fh.write(' {} {}\n'.format(n, m))

    # Print the first block.
    names = [name.replace(' ', '_') for name in names]
    name_length = max([len(name) for name in names]) + 2
    for name, seq in zip(names, seqs):
        fh.write(name + ' ' * (name_length - len(name)))
        fh.write(seq[:50])
        fh.write('\n')

    # Print the interleaved blocks.
    whitespace = ' ' * name_length
    i = 50
    while i <= m:
        fh.write('\n')
        for seq in seqs:
            fh.write(whitespace)
            fh.write(seq[i:i+50])
            fh.write('\n')
        i += 50

    fh.close()

def output_phylip_partitions(filename, partitions, alphabet):
    fh = open(filename, 'w')
    alphabet = alphabet.name.upper()
    template = alphabet+', Subset{} = {}\n'
    for i, part in enumerate(partitions):
        if part:
            # some bins may be empty
            l = template.format(i, str(part)[1:-1])
            fh.write(l)
    fh.close()

def output_nexus(filename, names, seqs, charsets, hist):
    # not yet implemented
    pass


def histogram(num_list, name_list):
    upper = float(max(num_list))
    pad = len(str(upper))
    parts = []
    for i in range(1,61):
        parts.append((upper/60)*i)

    for m, n in enumerate(num_list):
        pr = name_list[m] + "|"
        low = 0.0
        if n == 0:
            pr = pr + " "*61
        for p, hi  in enumerate(parts):
            if n > low and n <= hi:
                pr = pr + "="*(p+1) + (" "*(60 - p))
                break
            low = hi
        print "[" + pr + "|" + str(n) + " "*(pad-len(str(n))) + "]"


