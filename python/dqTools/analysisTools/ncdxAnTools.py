from __future__ import absolute_import

import os.path
import cPickle
import numpy as N
import ncdx

from .cache import cache
from .reg import register

@register
@cache('times.all')
def getAllTimes(base, pattern=None):
    dir = os.path.join(base,'out')
    series = ncdx.Series(dir)
    return N.array(series.keys())


class PseudoSwimmerArray(object):
    def __init__(self, **attrs):
        self.__dict__.update(attrs)
        return
    
    pass

@register
def getSwimmers(base, time):
    series = ncdx.Series(os.path.join(base, 'out'))
    swimmers = series[time].swimmers
    return PseudoSwimmerArray(r=swimmers.positions.data,
                              v=swimmers.velocity.data,
                              n=swimmers.orientation.data)

@register
def getParams(base):
    return cPickle.load(file(os.path.join(base, 'params')))

@register
def getFluid(base, time):
    series = ncdx.Series(os.path.join(base, 'out'))
    fluid = series[time].fluid
    return fluid.velocity.data
