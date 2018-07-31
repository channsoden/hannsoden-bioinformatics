#!/usr/bin/env python
import subprocess as sp

class ParserError(Exception):
    pass

def get_absolute_positions(fastafile, base=0):
    """Reads a fasta file and returns a dictionary containing the start positions
    of each sequence if they were concatenated together in order.
    Uses 0-based indexing by default."""
    from Bio import SeqIO
    fh = open(fastafile, 'r')
    records = SeqIO.parse(fh, 'fasta')
    d = {}
    genomelength = base
    for rec in records:
        d[rec.name] = genomelength
        genomelength += len(rec.seq)
    d['end'] = genomelength
    fh.close()
    return d

def get_scaffold_lengths(fastafile):
    """Opens and reads a fasta formated genome and returns a dictionary with sequence names as keys and the length of the sequence as values."""
    records = fasta_to_dict(fastafile)
    d = {key: len(rec) for key, rec in list(records.items())}
    return d

def seq_length_list(fastafile):
    """Returns a list of the lengths of the sequences in a fasta file in order."""
    data = open(fastafile, 'r').read()

    seqs = data[1:].split('>')
    pairs = [seq.split('\n', 1) for seq in seqs]
    lengths = [len(pair[1].replace('\n', '')) for pair in pairs]

    data.close()
    return lengths

def total_length(fastafile):
    """Returns the total length of sequences in a fasta file."""
    grep = sp.Popen(['grep', '-v', '">"', fastafile], stdout=sp.PIPE)
    wc = sp.Popen(['wc'], stdin=grep.stdout, stdout=sp.PIPE)
    out, err = wc.communicate()
    lines, words, characters = list(map(int, out.split()))
    return characters - lines

def fasta_to_dict(fastafile):
    """Very fast, but memory inefficient, function to return a dictionary of all sequences in the fasta file."""
    data = open(fastafile, 'r').read()

    genes = data[1:].split('>')
    genes = [gene.split('\n', 1) for gene in genes]
    try:
        genes = {gene[0] : gene[1].replace('\n', '') for gene in genes}
    except IndexError:
        raise ParserError('Empty or malformed sequences in {}'.format(fastafile))

    data.close()
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

    fh.close()
    out.close()

def fasta_to_phylip(fastafile):
    # Bio.AlignIO.write puts spaces in the sequence blocks that break some tools.
    records = fasta_to_dict(fastafile)
    seq_len = max([len(seq) for seq in list(records.values())])

    outfile = fastafile.split('/')[-1].rsplit('.', 1)[0] + '.phy'
    outfh = open(outfile, 'w')

    # Print the header
    outfh.write(' {} {}\n'.format(len(records), seq_len))
    
    # Print the first block.
    names = list(records.keys())
    name_length = max([len(name) for name in names]) + 2
    for name in names:
        outfh.write(name + ' ' * (name_length - len(name)))
        outfh.write(records[name][:50])
        outfh.write('\n')

    # Print the interleaved blocks.
    i = 50
    while i <= seq_len:
        outfh.write('\n')
        for name in names:
            outfh.write(' ' * name_length)
            outfh.write(records[name][i:i+50])
            outfh.write('\n')
        i += 50

    outfh.close()
    return outfile
