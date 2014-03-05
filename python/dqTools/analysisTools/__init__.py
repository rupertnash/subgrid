"""Some tools for analysing results of swimmers simulations.
"""
from __future__ import absolute_import

import os.path
from inspect import getmodule

from .cache import cache

import analysisTools.spAnTools
import analysisTools.ncdxAnTools
from .reg import __registry

__all__ = list(
    reduce(set.union,
           (_mfd.keys() for _mfd in __registry.itervalues()),
           set())
    )

def guessMod(base):
    ncOut = os.path.join(base, 'out')
    if os.path.exists(ncOut) and os.path.isdir(ncOut):
        return ncdxAnTools
    else:
        return spAnTools
    return

self = getmodule(guessMod)

for _fname in __all__:
    def definer(fname=_fname):
        def switcher(base, *args, **kwargs):
            mod = guessMod(base)
            func = __registry[mod.__name__][fname]
            return func(base, *args, **kwargs)
        return switcher
    
    setattr(self, _fname, definer())
    continue

__all__.append('cache')
