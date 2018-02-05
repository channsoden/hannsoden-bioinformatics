#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os

# Nonstandard modules

# My modules
from fasta_tools import total_length
from SLURM_tools import submit
from SLURM_tools import job_wait
import WGP_config as cfg

def partition(args, alignment):
    basedir = os.getcwd()
    try:
        os.mkdir('3_partitioning')
    except OSError:
        pass
    os.chdir('3_partitioning')

    phylip_file = args.output + '.phylip'
    partition_file = args.output + '.partitions.txt'

    # Rate partitioning usually takes under 0.5s per unique pattern
    # in the alignment. There is no fast way to know how many unique
    # patterns there will be in an alignment of a given size.
    pps = 0.05 # empirical guess of high patterns per site
    records = len(args.genomes)
    sites = total_length(alignment) / records
    estimated_patterns = pps * sites
    estimated_runtime = int(estimated_patterns * 0.5)
    minutes = (estimated_runtime / 60) +1

    command = 'tiger2 -in {} -a dna -out {} -f phylip -bt rota -b 4 -t 1'
    command = command.format(alignment, args.output)
    ID = submit(command,
                partition = cfg.SLURMpartition,
                account = cfg.SLURMaccount,
                qos = cfg.SLURMqos,
                time = str(minutes)
                job_name = 'rate_partitioning',
                cpus_per_task = cfg.SLURMcpus,
                mem_per_cpu = cfg.SLURMmem,
                modules = [cfg.python])
    job_wait(ID)
    outfile = 'rate_partitioning_'+str(ID)+'.out'
    errfile = 'rate_partitioning_'+str(ID)+'.err'

    os.chdir(basedir)
    cleanup([outfile, errfile])
    return basedir+'/3_partitioning/'+partition_file, basedir+'/3_partitioning/'+phylip_file

def cleanup(logs):
    try:
        os.mkdir('3_partitioning/logs')
    except OSError:
        pass
    for log in logs:
        os.rename('3_partitioning/'+log, '3_partitioning/logs/'+log)
