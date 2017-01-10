#!/usr/bin/env python

# submit("echo hello world", job_name='test', modules=['python', 'java', 'g++'])

import subprocess as sp
import time, sys

def add_option(options, option, value):
    if value:
        return '{} {}={}'.format(options, option, value)
    else:
        return options

def submit(job,
           partition = 'savio',
           account='co_rosalind',
           qos = 'rosalind_savio_normal',
           job_name = 'slurmjob',
           time = '48:0:0',
           error = '', 
           output = '', 
           nodes = '1',
           tasks_per_node = '1',
           cpus_per_task = '20',
           mem_per_cpu = '3000', 
           mail_type = 'FAIL',
           mail_user = 'no@nonexistent.nope',
           shell = 'bash',
           modules = [],
           verbose = True):

    if not error:
        error = job_name + '_%j.err'
    if not output:
        output = job_name + '_%j.out'

    options = {'--partition':partition, '--qos':qos, '--job-name':job_name,
               '--time':time, '--error':error, '--output':output,
               '--nodes':nodes, '--tasks-per-node':tasks_per_node,
               '--cpus-per-task':cpus_per_task, '--mem-per-cpu':mem_per_cpu,
               '--mail-type':mail_type, '--mail-user':mail_user, '--account':account}

    opt_string = ''
    for option, value in options.items():
        opt_string = add_option(opt_string, option, value)
    command = 'sbatch '+opt_string

    mod_string = ''
    for mod in modules:
        mod_string += 'module load '+mod+';\n'

    submission = sp.Popen(command, shell=True,
                          stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    job = '#!/bin/{}\n{}{}'.format(shell, mod_string, job)
    out, err = submission.communicate(input = job)
    if err:
        sys.exit(err)
    else:
        jobID = out.strip().split()[-1]
        if verbose:
            print 'Submitted: {}'.format(jobID)
        return jobID

def submit_script(job):
    submission = sp.Popen('sbatch ' + job, shell=True,
                          stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = submission.communicate()
    if err:
        sys.exit(err)
    else:
        jobID = out.strip().split()[-1]
        return jobID

def check_job(jobID):
    check = sp.Popen('sacct --format End -j {}'.format(jobID), shell=True,
                     stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = check.communicate()
    endtime = out.split('\n')[2].strip()
    finished = endtime != 'Unknown'
    return finished

def job_wait(jobID, period=60, verbose=True):
    unfinished = True
    while unfinished:
        # takes a second for the slurm scheduler to respond to sacct
        time.sleep(period)
        unfinished = not check_job(jobID)
    if verbose:
        print 'Completed: {}'.format(jobID)

