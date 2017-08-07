#!/usr/bin/env python
import sys
import fileinput

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

    while True:
        line = gffFH.readline()
        if not line.strip() or line.strip() == '##FASTA':
            break
        elif line[0] == '#':
            pass
        else:
            yield gff_feature(line)

