from __future__ import absolute_import

import os.path
import glob
import cPickle
import numpy as np

from .cache import cache

def getSwimmers(base, time):
    return cPickle.load(file(os.path.join(base,'swimmers.%.6d' % time)))

@cache('times.all')
def getAllTimes(base, pattern='swimmers.*'):
    gen_files = glob.iglob(os.path.join(base, pattern))
    gen_times = (int(f.rsplit('.', 1)[1]) for f in gen_files)
    times = sorted(gen_times)
    times = np.array(times, dtype=int)
    return times

def getParams(base):
    return cPickle.load(file(os.path.join(base, 'params')))
