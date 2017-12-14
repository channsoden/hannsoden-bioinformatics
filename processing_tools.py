#!/usr/bin/env python
"""
from processing_tools import mapPool

def f(arg1, ...argn):
    do_stuff
    return result

calls = [(f, (arg1, ...argn)) for item in list]

# Example for long/instense function:
results = mapPool(threads, calls, daemonic = True)

# Example for many short/small function calls:
results = mapPool(threads, calls, daemonic = True, chunksize = len(calls) / threads + 1)

# Example for functions that spend much of their time waiting for another process/event:
results = mapPool(threads, calls, daemonic = False)
"""

import sys, traceback
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

def mapPool(threads, args, daemonic=False, chunksize=None):
    # args = [(function, (arg1, arg2), {keyword1:arg1, keyword2:arg2}), nextjob, ...]
    if daemonic:
        p = Pool(threads)
    else:
        p = ThreadPool(processes=threads)

    output = p.map(star_func, args, chunksize=chunksize)
    p.close()
    p.join()
    return output

def star_func(arglist):
    func = arglist[0]
    args = arglist[1]
    try:
        kwargs = arglist[2]
        try:
            return func(*args, **kwargs)
        except:
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))
    except IndexError:
        try:
            return func(*args)
        except:
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))
