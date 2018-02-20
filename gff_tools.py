#!/usr/bin/env python
import sys, fileinput
from cStringIO import StringIO
import pandas as pd

class gff_feature:
    def __init__(self, gff_feature_line):
        fields = gff_feature_line.strip().strip(';').split('\t')
        self.seqname, self.source, self.feature = fields[:3]
        self.start, self.end = (int(x) for x in fields[3:5])
        try:
            self.score = int(fields[5])
        except ValueError:
            self.score = None
        self.strand, self.frame = fields[6:8]
        fields[8] = fields[8].split('#')[0].strip(';') # remove comments
        try:
            # atr=val format
            self.attributes = {atr:val for atr, val in [couple.split('=') for couple in fields[8].split(';')]}
        except ValueError:
            # atr "val" format, as produced by cufflinks
            self.attributes = {atr:val.strip('"') for atr, val in [couple.split() for couple in fields[8].split(';')]}

def parse_gff(gffFile):
    gffFH = open(gffFile, 'r')

    for line in gffFH:
        if not line.strip() or line.strip() == '##FASTA':
            break
        elif line[0] == '#':
            pass
        else:
            yield gff_feature(line)

def gff_table(gffFile):
    fh = open(gffFile, 'r')
    table_form = []
    for line in fh:
        if not line.strip() or line.strip() == '##FASTA':
            break
        elif line[0] == '#':
            pass
        else:
            line = line.split('#')[0].strip() # remove comments
            table_form.append(line)
    table_form = '\n'.join(table_form)
    table_form = StringIO(table_form)
    header = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(table_form, sep='\t', names=header)
    return df
