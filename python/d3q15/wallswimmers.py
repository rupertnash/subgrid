from __future__ import with_statement, absolute_import

from contextlib import contextmanager

from . import potentials
from .swimmers import SwimmerArray
from .CyDelta import delta
import numpy as N

class WallSwimmerArray(SwimmerArray):
    @staticmethod
    def calcSigma(L, **kwargs):
        # lengthscale of pot. should be a dist swimmer can move in a few
        # timesteps, choose 10.
        eta = L.tau_s / 3.
        return 10. *kwargs['P'] /(6.*N.pi * eta * kwargs['a'])
    
    def __init__(self, L, **kwargs):
        walls = kwargs.get('walls', True)
        if walls is False:
            raise ValueError('flag "walls" must be specified as True')
        
        forced = kwargs.get('forced', True)
        if walls is False:
            raise ValueError('flag "forced" must be specified as True')
        kwargs['forced'] = forced
        
        # Create the potential that keeps swimmers from the wall
        # (we need the length (sigma) later...)
        
        sigma = kwargs.get('sigma',
                           self.calcSigma(L,**kwargs))
        
        # depth should be that which will stop a swimmer at that distance
        epsilon = kwargs.get('epsilon',
                             sigma*kwargs['P'] / 24.)
        
        p = potentials.WallLJPotential(epsilon, sigma)

        # Now, set walls to be False & init from the superclass,
        # we set walls = False to skip the bounds checking.
        kwargs['walls'] = False
        SwimmerArray.__init__(self, L, **kwargs)
        # reset walls to True
        self.walls = True

        self.potentials.append(p)
        return
    __init__.__doc__ = SwimmerArray.__init__.__doc__

    @contextmanager
    def wallsFalse(self):
        self.walls = False
        try:
            yield
        finally:
            self.walls = True
    
    def applyMove(self, lattice, rdots):
        with self.wallsFalse():
            SwimmerArray.applyMove(self, lattice, rdots)
        return
    
    @classmethod
    def _interp(cls, lattice, rlist):
        v = N.zeros(rlist.shape)
        for delta3d, indexing in cls._delta_support(lattice, rlist):
            v += delta3d * lattice.u[indexing]
            continue
        
        return v
    
    @staticmethod
    def _delta_support(lattice,r):
        """Generator which iterates over the the support of Peskin's
        delta function, yielding the value of the delta function and
        """
        x = N.ceil(r[...,N.newaxis] + N.arange(-2,2).reshape(1,1,4))
        indices = x.astype(int)
        
        deltas = delta(x - r[...,N.newaxis])
        del x # don't need x anymore
        
        for d in range(3):
            # handle x & y as normal
            id = indices[:, d, :]
            n = lattice.size[d]
            
            where = N.where(id < 1)
            id[where] += n
            if d ==2:
                deltas[:,d,:][where] = 0.
                pass
            
            where = N.where(id>n)
            id[where] -= n
            if d==2:
                deltas[:,d,:][where] = 0.
                pass
            
            continue
        
        del where
        
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    # evaluate the 3d delta function 
                    delta3d = deltas[:,0:1,i] * deltas[:,1:2,j] * deltas[:,2:3,k]
                    yield (delta3d, (indices[:,0,i],
                                     indices[:,1,j],
                                     indices[:,2,k]))
                    continue
                continue
            continue
        return
    
    @classmethod 
    def _addForcesAtPositions(cls, lattice, force, r):
        for delta3d, indexing in cls._delta_support(lattice, r):
            lattice.force[indexing] += delta3d * force
            continue
        
    pass
