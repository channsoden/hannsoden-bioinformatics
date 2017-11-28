#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os, sys, pickle, shutil
import subprocess as sp

# Nonstandard modules

# My modules
from fasta_tools import fasta_to_dict
from SLURM_tools import submit
from SLURM_tools import job_wait

def partition(args, alignment):
    basedir = os.getcwd()
    try:
        os.mkdir('3_partitioning')
    except OSError:
        pass
    os.chdir('3_partitioning')

    try:
        os.mkdir(args.output)
    except OSError:
        pass
    
    phylip = fasta_to_phylip(alignment)

    # PartitionFinder must be run from a unique subdirectory
    # because the names of it's intermediate and output files
    # are hardcoded. This will allow for running multiple
    # instances of WGP from the same directory.
    phypath = basedir+'/3_partitioning/'+phylip
    link = basedir+'/3_partitioning/'+args.output+'/'+phylip
    if not (os.path.islink(link) and os.path.realpath(link) == phypath):
        try:
            os.remove(link)
        except OSError:
            pass
        os.symlink(phypath, link)
    os.chdir(args.output)
    configure_PF(phylip)
    os.chdir(basedir+'/3_partitioning')

    input_size = os.stat(phylip).st_size
    if input_size > 300 * 10 ** 6:
        ID = submit('{} {} {} 20'.format(sys.executable, __file__, basedir+'/'+args.output+'.args.pickle'),
                    partition = 'savio_bigmem',
                    account = 'co_rosalind',
                    qos = 'savio_lowprio',
                    time = '12:0:0',
                    job_name = 'partitionfinder',
                    cpus_per_task = 20,
                    mem_per_cpu = '25600',
                    modules = ['raxml/8.1.17'])
    elif input_size > 50 * 10 ** 6:
        ID = submit('{} {} {} 24'.format(sys.executable, __file__, basedir+'/'+args.output+'.args.pickle'),
                    partition = 'savio2_bigmem',
                    account = 'co_rosalind',
                    qos = 'savio_lowprio',
                    time = '12:0:0',
                    job_name = 'partitionfinder',
                    cpus_per_task = 24,
                    mem_per_cpu = '5300',
                    modules = ['raxml/8.1.17'])
    else:
        ID = submit('{} {} {} 20'.format(sys.executable, __file__, basedir+'/'+args.output+'.args.pickle'),
                    partition = 'savio',
                    account = 'co_rosalind',
                    qos = 'rosalind_savio_normal',
                    time = '12:0:0',
                    job_name = 'partitionfinder',
                    cpus_per_task = 20,
                    mem_per_cpu = '3000',
                    modules = ['raxml/8.1.17'])
    job_wait(ID)
    outfile = 'partitionfinder_'+str(ID)+'.out'
    errfile = 'partitionfinder_'+str(ID)+'.err'
    
    partition_file = get_scheme(args.output)

    os.chdir(basedir)
    cleanup(logs=[outfile, errfile], trashdir=basedir+'/3_partitioning/'+args.output)
    return basedir+'/3_partitioning/'+partition_file, phypath

def fasta_to_phylip(fastafile):
    # Had to write my own phylip writer because Bio.AlignIO.write puts spaces in the sequence blocks that break PartitionFinder.
    records = fasta_to_dict(fastafile)
    seq_len = max([len(seq) for seq in records.values()])

    outfile = fastafile.split('/')[-1].rsplit('.', 1)[0] + '.phy'
    outfh = open(outfile, 'w')

    # Print the header
    outfh.write(' {} {}\n'.format(len(records), seq_len))
    
    # Print the first block.
    names = records.keys()
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

def configure_PF(alignment_file, user_tree = '', branchlengths = 'linked', models='GTR+G', criteria = 'aicc', partition = '', search = 'kmeans'):
    # Create a partition_finder.cfg file
    cfg = open('partition_finder.cfg', 'w')

    cfg.write("# ALIGNMENT FILE #\n")
    cfg.write("alignment = {};\n".format(os.path.basename(alignment_file)))
    if user_tree:
        # Link the user tree into the working directory if necessary
        treebase = os.path.basename(user_tree)
        if not treebase in os.listdir('.'):
            os.symlink(user_tree, treebase)
        cfg.write("user_tree_topology = {};\n".format(treebase))
    cfg.write("\n")

    cfg.write("# BRANCHLENGTHS #\n")
    cfg.write("branchlengths = {};\n".format(branchlengths))
    cfg.write("\n")

    cfg.write("# MODELS OF EVOLUTION #\n")
    cfg.write("models = {};\n".format(models))
    cfg.write("model_selection = {};\n".format(criteria))
    cfg.write("\n")

    cfg.write("# DATA BLOCKS #\n")
    cfg.write("[data_blocks]\n")
    if partition:
        exit('configure_pf(): Configuring PF with a user defined partition is not yet implemented. Only kmeans algorithm is implimented at this point.')
    else:
        with open(alignment_file, 'r') as fh:
            genomesize = int(fh.readline().strip().split()[1])
        cfg.write("genome = 1-{};\n".format(genomesize))
    cfg.write("\n")

    cfg.write("# SCHEMES #\n")
    cfg.write("[schemes]\n")
    cfg.write("search = {};\n".format(search))
    cfg.write("\n")

    cfg.write("# user schemes (see manual)\n")
    
    cfg.close()

def partition_finder(args):
    # Run Partition Finder 2
    # Using more than one thread does not seem to make a difference, at least on my system with the current version.
    PFpath = '/global/scratch/hannsode/pkgs/partitionfinder/PartitionFinder.py'
    basedir = os.getcwd()
    os.chdir(args.output)
    command = '{0} {1} {2} --raxml'.format(sys.executable, PFpath, os.getcwd())
    partitioning = sp.Popen(command.split())
    os.chdir(basedir)
    return partitioning.wait()

def get_scheme(output):
    # Pulls the partitioning scheme suitable for RAxML/ExaML out of the results of a PartitionFinder analysis.
    with open(output+'/analysis/best_scheme.txt', 'r') as best_scheme:
        subsets = [line for line in best_scheme if line.startswith('DNA, Subset')]
    outfile = output+'.best_scheme.partition'
    with open(outfile, 'w') as fh:
        fh.writelines(subsets)
    return outfile

def cleanup(logs=[], trashdir=None):
    if logs and not os.path.isdir('3_partitioning/logs'):
        os.mkdir('3_partitioning/logs')
    [os.rename('3_partitioning/'+log, '3_partitioning/logs/'+log) for log in logs]
    if trashdir:
        shutil.rmtree(trashdir)

if __name__ == '__main__':
    argfile = sys.argv[1]
    with open(argfile, 'rb') as fh:
        args = pickle.load(fh)
    args.procs = int(sys.argv[2])

    partition_finder(args)
