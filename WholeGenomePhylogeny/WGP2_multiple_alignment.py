#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os, sys
import subprocess as sp
from multiprocessing import Pool

# Nonstandard modules
from Bio import SeqIO
import pandas as pd

# My modules
import fasta_tools
from processing_tools import mapPool
from SLURM_tools import submit
from SLURM_tools import job_wait

def multiple_alignment(args, fastas):
    basedir = os.getcwd()
    os.chdir('2_alignment')

    unaligned_fastas = [fasta for fasta in fastas if not os.path.isfile(trim_name(fasta))]
    if unaligned_fastas:
        chunk_size = (len(unaligned_fastas) / 4) + 1
        chunks = [fastas[i:i+chunk_size] for i in [n * chunk_size for n in range(4)]]
        # Run this script with list of fastas as args
        jobs = [(submit_alignment_batch, ['{} {} {}'.format(sys.executable, __file__, ' '.join(chunk))])
                for chunk in chunks]
        mapPool(4, jobs)

    aligned_trimmed = [trim_name(fasta) for fasta in fastas] # The output files from the aligment process.
    alignment = concatenate_fasta(aligned_trimmed, args.output)

    os.chdir(basedir)
    return basedir+'/2_alignment/'+alignment

def submit_alignment_batch(job):
    ID = submit(job,
                partition = 'savio',
                account = 'co_rosalind',
                qos = 'rosalind_savio_normal',
                time = '12:0:0',
                job_name = 'mafft',
                cpus_per_task = 20,
                mem_per_cpu = '3000')
                #modules = ['biopython/1.64-2.7.8'])
    job_wait(ID)

def align_trim(fasta):
    aligned = align_name(fasta)
    trimmed = trim_name(fasta)
    mafft = sp.Popen('mafft --globalpair --maxiterate 1000 --jtt 10 --nuc --inputorder {} 1> {} 2> /dev/null'.format(fasta, aligned), shell=True)
    mafft.wait()
    trim_gap_ends(aligned, trimmed)

def align_name(fasta):
    return os.path.splitext(fasta)[0] + '_aligned.fasta'

def trim_name(fasta):
    return os.path.splitext(fasta)[0] + '_aligned-trimmed.fasta'

def trim_gap_ends(infasta, outfasta):
    infh = open(infasta, 'r')
    outfh = open(outfasta, 'w')

    records = list(SeqIO.parse(infh, 'fasta'))
    boolArray = pd.DataFrame([list(rec.seq) for rec in records]) != '-'

    start = find_start(boolArray)
    end = find_end(boolArray)

    for record in records:
        record.seq = record.seq[start:end]

    SeqIO.write(records, outfh, 'fasta')
    infh.close()
    outfh.close()

def find_start(boolArray):
    upper = len(boolArray.columns)
    lower = 0
    while upper - lower > 1:
        i = (upper + lower) / 2
        if has_started(boolArray, i):
            upper = i
        else:
            lower = i
    return lower

def find_end(boolArray):
    upper = len(boolArray.columns)
    lower = 0
    while upper - lower > 1:
        i = (upper + lower) / 2
        if has_ended(boolArray, i):
            upper = i
        else:
            lower = i
    return upper

def has_started(boolArray, pos):
    return boolArray.iloc[:,:pos].any(axis=1).all()

def has_ended(boolArray, pos):
    return not boolArray.iloc[:,pos:].any(axis=1).all()

def concatenate_fasta(faAlignments, out_name):
    """Concatenates a list of fasta-formated alignment files."""
    for fa in faAlignments:
        try:
            seqs = fasta_tools.fasta_to_dict(fa)
            break
        except fasta_tools.ParserError:
            pass

    try:
        concatenated = {name.split()[-1].split(':')[0]: [] for name in seqs}
    except NameError:
        raise fasta_tools.ParserError('All alignments are empty or malformed.')
    
    bad_files = 0
    for fa in faAlignments:
        try:
            seqs = fasta_tools.fasta_to_dict(fa)
            [concatenated[name.split()[-1].split(':')[0]].append(seq)
             for name, seq in seqs.items()]
        except ParserError.ParserError:
            bad_files += 1    

    if bad_files:
        sys.stderr.write( 'Skipping {} empty or malformed alignment files\n'.format(bad_files) )
    
    # Abbreviate thxe sequence names and join all the sequences into a single long string.
    concatenated = {shorten_name(name): ''.join(seqs) for name, seqs in concatenated.items()}

    outfile = out_name+'.fasta'
    fh = open(outfile, 'w')
    for name, seq in concatenated.items():
        fh.write('>'+name+'\n')
        [fh.write(seq[i:i+80]+'\n') for i in range(0, len(seq), 80)]
    fh.close()
    return outfile

def shorten_name(sequence_name):
    return '_'.join(sequence_name.strip().split('/')[-1].split('.')[0].split('_')[:2])

def RAxML_valid(seqlist):
    ustates = set()
    for seq in seqlist:
        if len(seq) < 1:
            return False
        seq = seq.upper()
        ustates = ustates.union(set(seq))
        if len(ustates) >= 4:
            return True
    return False

if __name__ == '__main__':
    fastas = sys.argv[1:]
    calls = [(align_trim, [fasta]) for fasta in fastas]
    nones = mapPool(20, calls)
