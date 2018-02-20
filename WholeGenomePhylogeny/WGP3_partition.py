#!/usr/bin/env python
# Standard modules
import os, sys

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

    if (not os.path.isfile(phylip_file) or
        not os.path.isfile(partition_file)):
        submit_tiger2(args, alignment)

    os.chdir(basedir)
    return basedir+'/3_partitioning/'+partition_file, basedir+'/3_partitioning/'+phylip_file

def submit_tiger2(args, alignment):
    # Rate partitioning usually takes under 0.5s per unique pattern
    # in the alignment. There is no fast way to know how many unique
    # patterns there will be in an alignment of a given size.
    pps = 0.1 # empirical guess of high patterns per site # 0.05 for most, 0.1 for extremely diverse clade N. rajui
    records = len(args.genomes)
    sites = total_length(alignment) / records
    estimated_patterns = pps * sites
    estimated_runtime = int(estimated_patterns * 0.5)
    minutes = (estimated_runtime / 60) +1
    if minutes > (cfg.LARGEmaxtime * 60):
        warning = "Warning: estimated partition time ({}) is greater than maximum wallclock time ({}).\n"
        warning = warning.format(minutes, cfg.LARGEmaxtime * 60)
        sys.stderr.write(warning)
        minutes = cfg.LARGEmaxtime * 60

    command = '{} -in {} -a dna -out {} -f phylip -bt rota -b 4 -t 1'
    command = command.format(cfg.tiger, alignment, args.output)
    ID = submit(command,
                partition = cfg.LARGEpartition,
                account = cfg.LARGEaccount,
                qos = cfg.LARGEqos,
                time = str(minutes),
                job_name = 'rate_partitioning',
                cpus_per_task = cfg.LARGEcpus,
                mem_per_cpu = cfg.LARGEmem,
                modules = cfg.modules)
    job_wait(ID)

    outfile = 'rate_partitioning_'+str(ID)+'.out'
    errfile = 'rate_partitioning_'+str(ID)+'.err'
    cleanup([outfile, errfile])

def cleanup(logs):
    try:
        os.mkdir('logs')
    except OSError:
        pass
    for log in logs:
        os.rename(log, 'logs/'+log)
