import numpy as N

class Binner(object):
    """Class to perform averages of arbitary variables over some bin
    (i.e. range) of a coordinate.
    

    
    An example:
    
    r, u = calculateVelocitiesAtPositions()
    # r is nx X ny X nz X 3 array of positions
    # v is the same but for velocities at the corresponding position.
    
    # Distance from origin
    rNorm = N.sqrt(N.sum(r**2, axis=-1))
    
    # velocity magnitude
    uNorm = N.sqrt(N.sum(u**2, axis=-1))
    
    b = Binner(rNorm.ravel(), bins=100)
    
    uBin = b.mean(uNorm.ravel())
    
    g = Gnuplot.Gnuplot()
    g('set datafile missing "nan"')
    g('set logscale')
    
    g.plot(Gnuplot.Data(b.means,
                        b.mean(uNorm.ravel()),
                        b.map(N.std, rNorm.ravel()),
                        b.map(N.std, uNorm.ravel()),
                        **{'with': 'xyerrorlines'}))
    
    """

    def __init__(self, coord, **kwargs):
        """coord -- the coordinate on which we're constructing the bins
        
        The following arguments are as for numpy.histogram.
        
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width
            bins in the given range (10, by default). If `bins` is a sequence,
            it defines the bin edges, including the rightmost edge, allowing
            for non-uniform bin widths.
        range : (float, float), optional
            The lower and upper range of the bins.  If not provided, range
            is simply ``(a.min(), a.max())``.  Values outside the range are
            ignored. Note that with `new` set to False, values below
            the range are ignored, while those above the range are tallied
            in the rightmost bin.
        weights : array_like, optional
            An array of weights, of the same shape as `a`.  Each value in `a`
            only contributes its associated weight towards the bin count
            (instead of 1).  If `normed` is True, the weights are normalized,
            so that the integral of the density over the range remains 1.
        
        """
        self.coordShape = coord.shape
        
        self.counts, self.edges = N.histogram(coord, new=True, normed=False,
                                              **kwargs)
        self.centres = 0.5*(self.edges[:-1]+self.edges[1:])
        
        self.widths = self.edges[1:]-self.edges[:-1]
        self.nBins = len(self.counts)
        
        self.binmap = N.zeros(self.coordShape, dtype=N.int)
        for i in range(self.nBins):
            self.binmap[N.where(coord>self.edges[i])] = i
            continue
        
        self.masks = [N.ma.masked_not_equal(self.binmap, i, copy=False).mask
                      for i in range(self.nBins)]
        
        self.maskedCoords = [N.ma.array(coord, mask=msk)
                             for msk in self.masks]
        
        self.means = N.array([mc.mean() for mc in self.maskedCoords])
        
        return

    def map(self, func, data):
        """For each bin of the coordinate specified in the
        constructor, take an array of the corresponding elements from
        "data" and apply "func" to it, returning a numpy array of the
        result.
        
        func -- callable

        data -- numpy array of data, same shape as coord
        
        """
        return N.array([func(N.ma.array(data, mask=msk))
                        for msk in self.masks])
    
    def mean(self, data):
        """For each bin of the coordinate specified in the
        constructor, take an array of the corresponding elements from
        "data" and calculate it's mean, returning a numpy array of the
        results.
        
        data -- numpy array of data, same shape as coord
        
        """
        return self.map(N.mean, data)
    
    
if __name__ == '__main__':
    import cPickle
    import Gnuplot
    
    u = cPickle.load(file('/Disk/radio3data1/s0567077/pd_64/u.000122880'))
    uNorm = N.sqrt(N.sum(u**2, axis=-1))
    r = N.mgrid[(1.-u.shape[0])/2.:(1.+u.shape[0])/2.,
                (1.-u.shape[1])/2.:(1.+u.shape[1])/2.,
                (1.-u.shape[2])/2.:(1.+u.shape[2])/2.,].transpose((1,2,3,0))
    rNorm = N.sqrt(N.sum(r**2, axis=-1))

    b = Binner(rNorm.ravel(), bins=100, range=(0, u.shape[0]/2.))
    
    uBin = b.mean(uNorm.ravel())
    
    g = Gnuplot.Gnuplot()
    g('set datafile missing "nan"')
    g('set logscale')
    
    g.plot(Gnuplot.Data(b.means,
                        b.mean(uNorm.ravel()),
                        b.map(N.std, rNorm.ravel()),
                        b.map(N.std, uNorm.ravel()),
                        **{'with': 'xyerrorlines'}),
           '2e-4*x**-3')
    
    
