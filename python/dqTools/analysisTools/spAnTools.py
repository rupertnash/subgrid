from __future__ import absolute_import

import os.path
import glob
from swimPickle import c as cPickle
import numpy as N

from .cache import cache
from .reg import register

@register
def getSwimmers(base, time):
    return cPickle.load(file(os.path.join(base,'swimmers.%.9d' % time)))

@register
@cache('times.all')
def getAllTimes(base, pattern='swimmers.*'):
    gen_files = glob.iglob(os.path.join(base, pattern))
    gen_times = (int(f.rsplit('.', 1)[1]) for f in gen_files)
    times = sorted(gen_times)
    times = N.array(times, dtype=int)
    return times

@register
def getParams(base):
    return cPickle.load(file(os.path.join(base, 'params')))
