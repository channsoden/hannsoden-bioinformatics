#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os

# Nonstandard modules

# My modules
from SLURM_tools import submit
from SLURM_tools import job_wait

def partition(args, alignment):
    basedir = os.getcwd()
    try:
        os.mkdir('3_partitioning')
    except OSError:
        pass
    os.chdir('3_partitioning')

    input_size = os.stat(alignment).st_size / (10 ** 6)
    # make some guesses about runtime
    if input_size > 1000:
        time = '72:0:0'
    elif input_size > 300:
        time = '24:0:0'
    elif input_size > 150:
        time = '6:0:0'
    elif input_size > 50:
        time = '4:0:0'
    elif input_size > 10:
        time = '1:0:0'
    elif input_size > 1:
        time = '0:20:0'
    else:
        time = '0:1:0'
    
    command = 'tiger2 -in {} -a dna -out {} -f phylip -bt rota -b 4 -t 20'
    command = command.format(alignment, args.output)
    ID = submit(command,
                partition = 'savio',
                account = 'co_rosalind',
                qos = 'rosalind_savio_normal',
                time = time,
                job_name = 'rate_partitioning',
                cpus_per_task = 20,
                mem_per_cpu = '3000')
    job_wait(ID)
    outfile = 'rate_partitioning_'+str(ID)+'.out'
    errfile = 'rate_partitioning_'+str(ID)+'.err'

    phylip_file = args.output + '.phylip'
    partition_file = args.output + '.phartitions.txt'
    
    os.chdir(basedir)
    cleanup(logs=[outfile, errfile])
    return basedir+'/3_partitioning/'+partition_file, basedir+'/3_partitioning/'+phylip_file

def cleanup(logs=[]):
    if logs and not os.path.isdir('3_partitioning/logs'):
        os.mkdir('3_partitioning/logs')
    [os.rename('3_partitioning/'+log, '3_partitioning/logs/'+log) for log in logs]
