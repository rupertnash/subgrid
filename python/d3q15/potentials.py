"""Classes for managing a set of potentials. The potentials exert
forces which are a function of position.

"""

from __future__ import absolute_import
import numpy as N

class PotentialManager(list):
    """A container type for a list of potentials."""
    def __init__(self, seq=None):
        if seq is None:
            return list.__init__(self)
        else:
            seq = self._makeiterablecheck(seq)
            return list.__init__(self, seq)
        
    def __repr__(self):
        return self.__class__.__name__ + '(' + list.__repr__(self) + ')'

    
    def _check(self, item):
        if not isinstance(item, Potential):
            raise TypeError('All managed potentials must be an instance of (a subclass of) Potenial.Potential')

    def _makeiterablecheck(self, item):
        def check(item):
            for val in item:
                self._check(val)
                yield val
        return check(item)
            
    def __setitem__(self, key, value):
        if hasattr(value, '__iter__'):
            value = self._makeiterablecheck(value)
        else:
            self._checktype(value)
        
        return list.__setitem__(self, key, value)

    def __setslice__(self, i,j, seq):
        self.__setitem__(slice(i,j), seq)
        
    def append(self, value):
        self._check(value)
        return list.append(self, value)

    def extend(self, iterable):
        return list.extend(self, self._makeiterablecheck(iterable))
    
    def insert(self, index, item):
        self._check(item)
        return list.insert(self, index, item)

    def __add__(self, other):
        ans = self.__class__(self)
        ans.extend(other)
        return ans
    
    def __iadd__(self, other):
        self.extend(other)    

    def addForces(self, F, r, **kwargs):
        """Calculate the force due to all managed Potentials, at the
        points specified by 'r' and add it to the forces stored in
        'F'. The keyword arguments can be used.
        
        """
        for p in self:
            p.addForce(F, r, **kwargs)
            continue
        return

    pass

class Potential(object):
    """Base class for a potential."""
    def __call__(self, r, **kwargs):
        raise NotImplementedError('Potential must have form given')
    
    def force(self, r, **kwargs):
        raise NotImplementedError('Force must have form given')
    
    def addForce(self, F, r, **kwargs):
        F += self.force(r, **kwargs)
        return
    pass

class PeriodicBCPotential(Potential):
    """For cases with periodic boundary conditions. Subclasses are
    responsible for setting the origin attribute."""
    
    def minim(self, r, lattice):
        """Calculate the minimum image displacement from the origin."""
        dr = r - self.origin.reshape((r.ndim-1)*(1,)+(3,))
        for dim in range(3):
            ind = N.where(dr[:,dim] > +0.5 * lattice.size[dim])
            dr[ind, dim] -= lattice.size[dim]
            
            ind = N.where(dr[:,dim] < -0.5 * lattice.size[dim])
            dr[ind, dim] += lattice.size[dim]
            continue
        return dr
    pass

class HarmonicPotential(PeriodicBCPotential):
    """Spherical, harmonic potential."""
    def __init__(self, k, origin):
        """Potential is located at origin with spring constant k."""
        self.k = k
        self.origin = origin
        return
    
    def __call__(self, r, **kwargs):
        # return the potential at the positions r
        lattice = kwargs['lattice']
        dr = self.minim(r, lattice)
        return 0.5* self.k * N.sum(dr**2, axis=-1)
    
    def force(self, r, **kwargs):
        lattice = kwargs['lattice']
        dr = self.minim(r, lattice)
        return -self.k * dr
    pass

class HarmonicSlotPotential(PeriodicBCPotential):
    
    def __init__(self, k, origin, normal):
        self.k = k
        self.origin = origin
        self.normal = normal / N.sqrt(N.sum(normal**2))
        
        return

    def __call__(self, r, **kwargs):
        lattice = kwargs['lattice']
        dr = self.minim(r, lattice)
        rnotnsq= N.sum(dr * self.normal.reshape((r.ndim-1)*(1,)+(3,)), axis=-1)
        
        return 0.5 * self.k * rdotnsq

    def force(self, r, **kwargs):
        lattice = kwargs['lattice']
        dr = self.minim(r, lattice)
        
        rnotnsq= N.sum(dr * self.normal.reshape((r.ndim-1)*(1,)+(3,)), axis=-1)
        return -self.k * self.normal.reshape((r.ndim-1)*(1,)+(3,)) * \
               N.sqrt(rdotnsq)
    pass

class WallLJPotential(Potential):
    """A Lennard-Jones type potential for keeping swimmers away from walls.
    This is rule B from my thesis."""
    
    def __init__(self, epsilon, sigma, walldim=2):
        self.epsilon = epsilon
        self.sigma = sigma
        self.rStar =  sigma * 2**(1./6.)
        self.walldim = walldim
        return
    
    def __call__(self, r, **kwargs):
        lattice = kwargs['lattice']
        z = r[..., self.dim] - 0.5
        
        U = N.zeros(z.shape)
        
        lowers = N.where(z < self.rStar)
        sigma_z_6 = (self.sigma/z[lowers])**6
        U[lowers] = (4*self.epsilon)*(sigma_z_6**2 - sigma_z_6) + self.epsilon
        del lowers
        
        uppers = N.where(z > lattice.size[self.dim]-self.rStar)
        sigma_z_6 = (self.sigma/(lattice.size[self.dim] - z[uppers]))**6
        U[uppers] = (4*self.epsilon)*(sigma_z_6**2 - sigma_z_6) + self.epsilon
        return U

    def force(self, r, **kwargs):
        lattice = kwargs['lattice']
        z = r[..., self.walldim] - 0.5
        
        F = N.zeros(r.shape)
        Fz = F[...,self.walldim]
        
        lowers = N.where(z < self.rStar)
        dz = z[lowers]
        sigma_z_6 = (self.sigma/dz)**6
        Fz[lowers] = (24*self.epsilon)*(2*sigma_z_6**2 - sigma_z_6) / \
                     dz
        del lowers
        
        uppers = N.where(z > lattice.size[self.walldim]-self.rStar)
        dz = z[uppers] - lattice.size[self.walldim]
        sigma_z_6 = (self.sigma/dz)**6
        Fz[uppers] = (24*self.epsilon)*(2*sigma_z_6**2 - sigma_z_6) / dz
        return F
