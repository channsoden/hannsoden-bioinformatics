#!/usr/bin/env python

# from SLURM_tools import submit, job_wait
# JID = submit("echo hello world", job_name='test', modules=['python', 'java', 'g++'])
# job_wait(JID, period = 15, verbose = False)

import subprocess as sp
import time, sys, tempfile

def add_option(options, option, value):
    if value:
        return '{} {}={}'.format(options, option, value)
    else:
        return options

def submit(job,
           pool = False,
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
           mail_user = 'channsoden@berkeley.edu',
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

    # pool allows submitting a series of Python function calls to be run in parrallel
    # job = [(function, (arg1, arg2), {keyword1:arg1, keyword2:arg2}), nextjob, ...]
    if pool:
        this_python = sys.executable
        script = ['from processing_tools import mapPool\n']

        # Need to convert job into string of python code, using f.__name__ without quotes
        try:
            module_functions = [(sys.modules[f.__module__].__file__.rsplit('.', 1)[0].split('/')[-1],
                                 f.__name__) for f, args in job]
            job = '['+', '.join(['({}, {})'.format(f.__name__, args) for f, args in job])+']'
        except ValueError:
            try:
                module_functions = [(sys.modules[f.__module__].__file__.rsplit('.', 1)[0].split('/')[-1],
                                     f.__name__) for f, args, kwargs in job]
                job = '['+', '.join(['({}, {}, {})'.format(f.__name__, args, kwargs) for f, args, kwargs in job])+']'
            except ValueError:
                raise ValueError('job list in pool should be all 2-uples (function, (args)) or all 3-uples (function, (args), {kwargs})')

        module_functions = set(module_functions)
        for module, func in module_functions:
            script.append('from {} import {}\n'.format(module, func))
        script.append('mapPool({}, {})'.format(cpus_per_task, job))

        temp_script = tempfile.NamedTemporaryFile(prefix = job_name+'_pool', dir='.', delete=False)
        temp_script.writelines(script)
        temp_script.close()
        job = '{} {}'.format(this_python, temp_script.name)

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
    if endtime:
        finished = endtime != 'Unknown'
    else:
        finished = False
    return finished

def job_wait(jobID, period=60, verbose=True):
    unfinished = True
    # takes a second for the slurm scheduler to respond to sacct
    time.sleep(1)
    while unfinished:
        unfinished = not check_job(jobID)
        time.sleep(period)
    if verbose:
        print 'Completed: {}'.format(jobID)


