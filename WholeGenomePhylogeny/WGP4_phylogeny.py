#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os, random, sys
import subprocess as sp

# Nonstandard modules

# My modules
from SLURM_tools import submit
from SLURM_tools import job_wait

def phylogeny(args, partition_file, phylip):
    random.seed(args.seed)

    basedir = os.getcwd()
    try:
        os.mkdir('4_phylogeny')
    except OSError:
        pass
    os.chdir('4_phylogeny')
    
    # Generate a parsimony starting tree using RAxML
    ST = starting_tree(phylip)

    # Parse the original data into an ExaML binary file
    partitioned = parse_examl(phylip, partition_file)

    # Perform an ExaML analysis on the original data to generate a tree.
    tree = examl(partitioned, ST, args.output)

    os.chdir(basedir)
    return basedir+'/4_phylogeny/'+tree
    
def starting_tree(alignment_file):
    # Compute parsimony starting trees for each replicate in same manner as intial inference.
    outfile = alignment_file.split('/')[-1] + '.ST'
    randomseed = random.randint(1, 999999)
    command = 'raxmlHPC-SSE3 -y -m GTRGAMMA -p {0} -s {1} -n {2}'.format(randomseed, alignment_file, outfile)
    parsimony = sp.Popen(command.split(), stdout=open('startingtree.out', 'a'), stderr=open('startingtree.err', 'a'))
    parsimony.wait()
    return 'RAxML_parsimonyTree.' + outfile

def parse_examl(alignment_file, partition_file):
    # Make an ExaML binary file from an alignment and partition
    outfile = alignment_file.split('/')[-1] + '.partitioned'
    outlog = 'parse_{}.out'.format(alignment_file.split('/')[-1].rsplit('.', 1)[0])
    parser_path = '/global/scratch/hannsode/pkgs/ExaML-master/parser/parse-examl'
    command = '{} -s {} -m DNA -q {} -n {}'.format(parser_path, alignment_file, partition_file, outfile)
    parsing = sp.Popen(command.split(), stdout=open(outlog, 'a'))
    parsing.wait()
    return outfile + '.binary'

def examl(binary_alignment, startingtree, outprefix):
    randomseed = random.randint(1, 999999)
    examl_path = '/global/scratch/hannsode/pkgs/ExaML-master/examl/examl-OMP-AVX'
    command = '{} -s {} -n {} -m GAMMA -t {} -p {}'.format(examl_path, binary_alignment, outprefix, startingtree, randomseed)
    outlog = 'examl_{}.out'.format(outprefix)
    errlog = 'examl_{}.err'.format(outprefix)
    sp.Popen(command, shell=True, stdout=open(outlog, 'a'), stderr=open(errlog, 'a')).wait()
    return 'ExaML_result.' + outprefix
