#!/usr/bin/env python
# Standard modules
import os, random, sys
import subprocess as sp

# Nonstandard modules

# My modules
import WGP_config as cfg
from SLURM_tools import submit, job_wait

def phylogeny(args, partition_file, phylip):
    basedir = os.getcwd()
    try:
        os.mkdir('4_phylogeny')
    except OSError:
        pass
    os.chdir('4_phylogeny')

    os.symlink(partition_file, partition_file.split('/')[-1])
    os.symlink(phylip, phylip.split('/')[-1])
    partition_file = partition_file.split('/')[-1]
    phylip = phylip.split('/')[-1]

    tree = partition_file + '.treefile'

    if not os.path.isfile(tree):
        job = '{} -nt AUTO -s {} -m MFP -alrt 1000 -bb 1000 -spp {}'
        job = job.format(cfg.iqtree, phylip, partition_file)
        ID = submit(job,
                    partition = cfg.SLURMpartition,
                    account = cfg.SLURMaccount,
                    qos = cfg.SLURMqos,
                    time = '24:0:0',
                    job_name = 'IQTree',
                    cpus_per_task = cfg.SLURMcpus,
                    mem_per_cpu = cfg.SLURMmem,
                    modules = cfg.modules)
        job_wait(ID)

        err = 'IQTree_{}.err'.format(ID)
        out = 'IQTree_{}.out'.format(ID)
        cleanup(logs=[out, err])

    os.chdir(basedir)
    return '4_phylogeny/'+tree

def cleanup(logs=[], trash=[]):
    try:
        os.mkdir('logs')
    except OSError:
        pass
    for log in logs:
        os.rename(log, 'logs/'+log)
    for f in trash:
        os.remove(f)
