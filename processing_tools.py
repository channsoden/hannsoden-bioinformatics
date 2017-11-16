#!/usr/bin/env python
"""
from processing_tools import mapPool

def f(arg1, ...argn):
    do_stuff
    return result

calls = [(f, (arg1, ...argn)) for item in list]
results = mapPool(threads, calls)
"""
import sys, traceback
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

def mapPool(threads, args, daemonic=False):
    # args = [(function, (arg1, arg2), {keyword1:arg1, keyword2:arg2}), nextjob, ...]
    if daemonic:
        p = Pool(threads)
    else:
        p = ThreadPool(processes=threads)

    output = p.map(star_func, args)
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
    
