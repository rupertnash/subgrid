"""Quick creation of histograms for plotting in Gnuplot.
"""
import numpy as N
import Gnuplot

def autohist(data, histopts={}, plotopts={}):
    dens, edges = N.histogram(data, **histopts)
    x = 0.5 * (edges[:-1] + edges[1:])
    return Gnuplot.Data(x, dens, **plotopts)
