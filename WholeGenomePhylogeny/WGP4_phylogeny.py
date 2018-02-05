#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os, random, sys
import subprocess as sp

# Nonstandard modules

# My modules
import WGP_config as cfg
from SLURM_tools import submit
from SLURM_tools import job_wait

def phylogeny(args, partition_file, phylip):
    picklefile = args.output+'.args.pickle'
    job = '{} {} {} {}'.format(sys.executable, __file__, picklefile, partition_file, phylip)
    submit_phylogeny(job)

    tree = os.getcwd() + '/4_phylogeny/ExaML_result.' + args.output
    return tree

def submit_phylogeny(job):
    ID = submit(job,
                partition = cfg.SLURMpartition,
                account = cfg.SLURMaccount,
                qos = cfg.SLURMqos,
                time = '12:0:0',
                job_name = 'ExaML',
                cpus_per_task = cfg.SLURMcpus,
                mem_per_cpu = cfg.SLURMmem,
                modules = [cfg.python, cfg.gcc, cfg.mpi])
    job_wait(ID)
    return ID

def main(args, partition_file, phylip)
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
    command = '{} -y -m GTRGAMMA -p {0} -s {1} -n {2}'.format(cfg.raxml, randomseed, alignment_file, outfile)
    parsimony = sp.Popen(command.split(), stdout=open('startingtree.out', 'a'), stderr=open('startingtree.err', 'a'))
    parsimony.wait()
    return 'RAxML_parsimonyTree.' + outfile

def parse_examl(alignment_file, partition_file):
    # Make an ExaML binary file from an alignment and partition
    outfile = alignment_file.split('/')[-1] + '.partitioned'
    outlog = 'parse_{}.out'.format(alignment_file.split('/')[-1].rsplit('.', 1)[0])
    command = '{} -s {} -m DNA -q {} -n {}'.format(cfg.parse_examl, alignment_file, partition_file, outfile)
    parsing = sp.Popen(command.split(), stdout=open(outlog, 'a'))
    parsing.wait()
    return outfile + '.binary'

def examl(binary_alignment, startingtree, outprefix):
    randomseed = random.randint(1, 999999)
    command = '{} {} -n {} -m GAMMA -t {} -p {}'.format(cfg.examl, binary_alignment, outprefix, startingtree, randomseed)
    outlog = 'examl_{}.out'.format(outprefix)
    errlog = 'examl_{}.err'.format(outprefix)
    outfh = open(outlog, 'w')
    errfh = open(errlog, 'w')
    sp.Popen(command, shell=True, stdout=outfh, stderr=errfh).wait()
    outfh.close()
    errfh.close()
    result = 'ExaML_result.' + outprefix

    partlog = 'parse_{}.out'.format(binary_alignment.rsplit('.', 3)[0])
    cleanup(logs=[outlog, errlog, partlog],
            trash=[f for f in os.listdir('.')
                   if (outprefix in f and
                       f  not in [result, outlog, errlog, partlog])])
    return result

def cleanup(logs=[], trash=[]):
    try:
        os.mkdir('logs')
    except OSError:
        pass
    for log in logs:
        os.rename(log, 'logs/'+log)
    for f in trash:
        os.remove(f)

if __name__ == '__main__':
    with open(sys.argv[1], 'rb') as fh:
        args = pickle.load(fh)

    part_file = sys.argv[2]
    phylip = sys.argv[3]
    tree = main(args, partition_file, phylip)
