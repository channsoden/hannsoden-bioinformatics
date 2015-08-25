#!/usr/bin/env python

def get_absolute_positions(fastafile):
    """Opens and reads a fasta formated genome and returns a dictionary suitable for use in parse_mummer_coords()."""
    from Bio import SeqIO
    records = SeqIO.parse(open(fastafile, 'r'), 'fasta')
    d = {}
    genomelength = 0
    for rec in records:
        d[rec.name] = genomelength
        genomelength += len(rec.seq)
    d['end'] = genomelength
    return d

def get_scaffold_lengths(fastafile):
    """Opens and reads a fasta formated genome and returns a dictionary with sequence names as keys and the length of the sequence as values."""
    from Bio import SeqIO
    records = SeqIO.parse(open(fastafile, 'r'), 'fasta')
    d = {rec.name: len(rec.seq) for rec in records}
    return d

def fasta_to_dict(fastafile):
    fh = open(fastafile, 'r')
    d = {}
    for line in fh:
        line = line.strip()
        if line[0] == '>':
            seqname = line[1:]
            d[seqname] = ''
        else:
            d[seqname] += line
    return d
