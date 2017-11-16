#!/usr/bin/env python
import subprocess as sp

def get_absolute_positions(fastafile, base=0):
    """Reads a fasta file and returns a dictionary containing the start positions
    of each sequence if they were concatenated together in order.
    Uses 0-based indexing by default."""
    from Bio import SeqIO
    records = SeqIO.parse(open(fastafile, 'r'), 'fasta')
    d = {}
    genomelength = base
    for rec in records:
        d[rec.name] = genomelength
        genomelength += len(rec.seq)
    d['end'] = genomelength
    return d

def get_scaffold_lengths(fastafile):
    """Opens and reads a fasta formated genome and returns a dictionary with sequence names as keys and the length of the sequence as values."""
    records = fasta_to_dict(fastafile)
    d = {key: len(rec) for key, rec in records.items()}
    return d

def seq_length_list(fastafile):
    """Returns a list of the lengths of the sequences in a fasta file in order."""
    data = open(fastafile, 'r').read()

    seqs = data[1:].split('>')
    pairs = [seq.split('\n', 1) for seq in seqs]
    lengths = [len(pair[1].replace('\n', '')) for pair in pairs]

    return lengths

def total_length(fastafile):
    """Returns the total length of sequences in a fasta file."""
    grep = sp.Popen(['grep', '-v', '">"', fastafile], stdout=sp.PIPE)
    wc = sp.Popen(['wc'], stdin=grep.stdout, stdout=sp.PIPE)
    out, err = wc.communicate()
    lines, words, characters = map(int, out.split())
    return characters - lines

def fasta_to_dict(fastafile):
    """Very fast, but memory inefficient, function to return a dictionary of all sequences in the fasta file."""
    data = open(fastafile, 'r').read()

    genes = data[1:].split('>')
    genes = [gene.split('\n', 1) for gene in genes]
    try:
        genes = {gene[0] : gene[1].replace('\n', '') for gene in genes}
    except IndexError:
        raise ValueError('Empty or malformed sequences in {}'.format(fastafile))

    return genes

def wrap_fasta(fasta_file, out_file = '', N = 80):
    """Wraps long sequences and writes the reformatted file back to [fasta_file].wrapped."""
    if not out_file:
        out_file = fasta_file + '.wrapped'

    fh = open(fasta_file, 'r')
    out = open(out_file, 'w')

    for line in fh:
        if line.startswith('>'):
            out.write(line)
        else:
            line = line.strip()
            for i in range(0, len(line), N):
                out.write(line[i: i+N])
                out.write('\n')
