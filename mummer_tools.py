#!/usr/bin/env python

class MUMmer_Coord:
    def __init__(self, refstart, refend, qstart, qend, refalen, qalen, ID, reflen, qlen, refcov, qcov, refscaf, qscaf, idx):
        self.refstart = refstart
        self.refend = refend
        self.qstart = qstart
        self.qend = qend
        self.refalen = refalen
        self.qalen = qalen
        self.ID = ID
        self.reflen = reflen
        self.qlen = qlen
        self.refcov = refcov
        self.qcov = qcov
        self.refscaf = refscaf
        self.qscaf = qscaf
        self.idx = idx

        if (refstart < refend and qstart < qend) or (refstart > refend and qstart > qend):
            self.orientation = True
        else:
            self.orientation = False

def parse_mummer_coords(showcoords_output, ref_scaffold_sizes,
                query_scaffold_sizes, header=True):
    """
    showcoords_output should be filehandle containting the output from
    show-coords with the -c and -l output options
    only. ref_scaffold_sizes and query_scaffold_sizes are dictionaries
    containing the names of all the scaffolds in the reference and
    query genomes as keys, and the absolute start position of the
    scaffolds as values. For example, scaffold1 starts at 0 and is
    1000bp long, then scaffold2 starts at 1000, etc. These
    dictionaries should also contain one key "end" which gives the
    total size of the genome.
    
    Creates a generator that yields one MUMer_Coord object for every
    alignment line in showcoords_output.
    """

    import numpy as np

    if header:
        while not set(showcoords_output.readline().strip()) == set('='): pass

    coordinate_array = []
    line = showcoords_output.readline()

    idx = 0
    while line:
        line = line.strip().replace('|','', 6).split()
        refstart, refend, qstart, qend, refalen, qalen = list(map(int, line[:6]))
        ID = float(line[6])
        reflen, qlen = list(map(int, line[7:9]))
        refcov, qcov = list(map(float, line[9:11]))
        refscaf, qscaf = line[11:]
        idx += 1
        
        # Convert from the relative positions (along scaffolds) to absolute
        # positions in the genome.
        refstart += ref_scaffold_sizes[refscaf]
        refend += ref_scaffold_sizes[refscaf]
        qstart += query_scaffold_sizes[qscaf]
        qend += query_scaffold_sizes[qscaf]

        coord = MUMmer_Coord(refstart, refend, qstart, qend, refalen,
                             qalen, ID, reflen, qlen, refcov, qcov,
                             refscaf, qscaf, idx)
        yield coord

        line = showcoords_output.readline()

class alignment:
    def __init__(self, rSeq, qSeq, rSeqLen, qSeqLen, rStart, rEnd, qStart, qEnd, errors, simErrors, stops, delta):
        self.rSeq = rSeq
        self.qSeq = qSeq
        self.rSeqLen = rSeqLen
        self.qSeqLen = qSeqLen
        self.rStart = rStart
        self.rEnd = rEnd
        self.qStart = qStart
        self.qEnd = qEnd
        self.errors = errors
        self.simErrors = simErrors
        self.stops = stops
        self.delta = delta

        self.rLen = rEnd - rStart
        self.qLen = abs(qEnd - qStart)

        self.deletions = len([digit for digit in delta if digit < 0])
        self.insertions = len([digit for digit in delta if digit > 0])
        self.aLen = self.rLen + self.deletions

        if qStart < qEnd: # Should rStart < rEnd, since it's the reference.
            # Same orientation is defined as True
            self.orientation = True
        else:
            self.orientation = False

def parse_delta(deltafile):
    fh = open(deltafile, 'r')
    fh.readline() # Input files
    fh.readline() # Type (NUCMER/PROMER)
    
    line = fh.readline()
    while line:
        if line[0] == '>':
            # This is a header
            spLine = line[1:].strip().split()
            if len(spLine) == 4:
                rSeq, qSeq = spLine[:2]
                rSeqLen, qSeqLen = list(map(int, spLine[2:]))
            else:
                exit('Invalid delta format. There may be white space in sequence names. Check your fasta files to make sure there is no white space used as delimeters.')
        else:
            # This is the start of an alignment
            rStart, rEnd, qStart, qEnd, errors, simErrors, stops = list(map(int, line.strip().split()))

            delta = []
            while True:
                digit = int(fh.readline().strip())
                if digit:
                    delta.append(digit)
                else:
                    break

            yield alignment(rSeq, qSeq, rSeqLen, qSeqLen, rStart, rEnd, qStart, qEnd, errors, simErrors, stops, delta)

        line = fh.readline()

def sorted_count(query, li):
    count = 0
    if li[0] == query:
        try: count += 1 + sorted_count(query, li[1:])
        except IndexError: count = 1
    return count

def map_coverage(coords, genome_length, perspective=''):
    """
    coords should be an iterable object that yields or contains
    MUMer_Coord objects, such as a mummer_tools.parse_mummer_coords
    generator. perspective should be either 'query' or 'reference',
    and genome_length should be the length of the query or reference
    genome, respectively.

    I want to annotate the genome of interest with the coverage of
    alignments to the opposite genome.
    [[start, end, coverage], [start, end, coverage]. . .]
    coverage is the number of times the genome of interest aligns to
    the other genome between start and end.
    """
    starts = []
    ends = []

    if perspective == 'query':
        for coord in coords:
            start = min(coord.qstart, coord.qend)
            end = max(coord.qstart, coord.qend)
            starts.append(start)
            ends.append(end)
    elif perspective == 'reference':
        for coord in coords:
            start = min(coord.refstart, coord.refend)
            end = max(coord.refstart, coord.refend)
            starts.append(start)
            ends.append(end)
    else:
        raise ValueError("Valid values for keyword 'perspective' in map_coverage() are 'query' or 'reference' only.")

    starts.sort()
    ends.sort()

    pos = 1
    coverage = 0
    cov_map = []
    segment_start = 1
    write = False
    delta_coverage = 0
    while pos < genome_length:
        if starts and pos == starts[0]:
            # Encountered the start of a segment.
            start_count = sorted_count(pos, starts)
            delta_coverage += start_count
            for x in range(start_count): starts.remove(pos)
            write = True
        if ends and pos == ends[0]:
            # Encountered the end of a segment.
            end_count = sorted_count(pos, ends)
            delta_coverage -= end_count
            for x in range(end_count): ends.remove(pos)
            write = True
        if write:
            cov_map.append([segment_start, pos-1, coverage])
            segment_start = pos
            coverage += delta_coverage
            delta_coverage = 0
            write = False
        pos += 1
    
    return cov_map
