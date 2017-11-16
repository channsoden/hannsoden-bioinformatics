#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os, shutil
import subprocess as sp

def bundle(glob, directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
    command = 'mv {} {}'.format(glob, directory)
    sp.Popen(command, shell=True).wait()

def cleanup(args):
    # Remove pickles
    [os.remove(f) for f in os.listdir('.') if f.endswith('.pickle')]
    [os.remove('1_orthology/'+f) for f in os.listdir('1_orthology') if f.endswith('.pickle')]

    # Bundle logs
    bundle('py_ortho*.???', '1_orthology/logs')
    bundle('1_orthology/mummer*.???', '1_orthology/logs')
    bundle('2_alignment/mafft*.???', '2_alignment/logs')
    bundle('3_partitioning/partitionfinder_*.???', '3_partitioning/logs')
    bundle('4_phylogeny/*.err 4_partitioning/*.out', '4_phylogeny/logs')

    # Remove alignments
    [os.remove('2_alignment/'+f) for f in os.listdir('2_alignment') if f.endswith('.fa') or f.endswith('.fasta')]

    # Remove partitionfinder files
    shutil.rmtree('3_partitioning/analysis/')
    [os.remove('3_partitioning/'+f) for f in os.listdir('3_partitioning') if not os.path.isdir('3_partitioning/'+f)]

    # Remove examl files
    [os.remove('4_phylogeny/'+f) for f in os.listdir('4_phylogeny') if
     (not os.path.isdir('4_phylogeny/'+f) and not f.startswith('ExaML_result.'))]

    
