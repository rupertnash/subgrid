"""Main swimmers class.

See my thesis for general discussion."""

from __future__ import absolute_import

import d3q15
from . import potentials, utils
import numpy as N

from dqTools import randomUnitVector

ceil = N.ceil

def norm(v):
    """Return the norm of a vector. Quick and dirty, don't use unless
    you're part of the package."""
    return N.sqrt(N.sum(v**2,axis=-1))


class SwimmerArray(d3q15.XArray):
    def __init__(self, L, **kwargs):
        """Create an array of Swimmer objects.

        Arguments:
          L: the Lattice object on which the Swimmers are moving

        Keyword arguments:
          a: the hydrodynamic radius of the particles
          l: the length of the swimmer
          r_list: the initial positions of the particles
          n_list: the initial orientations of the particles
          P: the magnitude of the force that is propelling the particles
          hydroRadius: the lattice hydrodynamic radius; see fallers.calibrate

        Optional keyword arguments:
          randSeed: seed for the PRNG. Zero => init from /dev/urandom
          walls: flag indicating whether the Lattice has walls
          tumbles: flag indicating whether the particles tumble
          tumbleProb: required if tumble given
          sinks: flag indicating an external force
          F: required if sinks given; the external force
          tracks: flag indicating whether to track the unwrapped position

        """
        try:
            a = kwargs['a']
            l = kwargs['l']
            P = kwargs['P']
            r_list = kwargs['r_list']
            n_list = kwargs['n_list']
            self.hydroRadius = kwargs['hydroRadius']
            assert self.hydroRadius > kwargs['a']
        except KeyError, err:
            raise TypeError("Must specify %s" % err)

        try:
            r_list = N.array(r_list)
            n_list = N.array(n_list)
        except ValueError, err:
            raise TypeError('%s must be or be convertable to a numpy array of floats' % err)
        
        try:
            randSeed = kwargs['randSeed']
        except KeyError:
            randSeed = 0

        # Set default options
        self.__setattr__('swims', True)
        self.__setattr__('advects', True)
        self.__setattr__('rotates', True)
        self.__setattr__('tumbles', True)
        self.__setattr__('tracks', True)
        self.__setattr__('walls', False)
        self.__setattr__('sinks', False)
        self.__setattr__('forced', False)

        for flag in ['swims', 'advects', 'rotates', 'tumbles', 'tracks', 'walls', 'sinks', 'forced']:
            try:
                if kwargs[flag]:
                    self.__setattr__(flag, True)
                else:
                    self.__setattr__(flag, False)
            except KeyError:
                continue

        if self.tumbles:
            try: 
                self.tumbleProb = kwargs['tumbleProb']
            except KeyError, err:
                raise TypeError("Must specify %s" % err)

        if self.sinks:
            try:
                F = kwargs['F']
            except KeyError, err:
                raise TypeError("Must specify %s" % err)

        self.eta = L.tau_s / 3.
        self.num = len(r_list)
 
        # set up fields in data array
        # we only create a table if one doesn't exist already
        # (we might be called by a subclass that created some already)
        if not hasattr(self, 'table'):
            self.table = {}
        self.tableAppend('r',3)
        if self.tracks:
            self.tableAppend('s', 3)
        self.tableAppend('v',3)
        self.tableAppend('n',3)
        self.tableAppend('P',1)
        self.tableAppend('a',1)
        self.tableAppend('l',1)
        if self.tumbles:
            self.tableAppend('t',1)
        if self.sinks:
            self.tableAppend('F',3)
        if self.forced:
            self.tableAppend('G',3)
        
        d3q15.XArray.__init__(self, L)
        
        try:
            # check positions in range
            assert N.all(r_list > 0.5)
            assert N.all(r_list[:, 0] <= L.nx+0.5)
            assert N.all(r_list[:, 1] <= L.ny+0.5)
            assert N.all(r_list[:, 2] <= L.nz+0.5)
            if self.walls:
                assert N.all(r_list[:,2] > 2.)
                assert N.all(r_list[:,2] < L.nz-1.)
        except AssertionError:
            raise ValueError("Coordinate(s) out of range")

        self.r = r_list
        if self.tracks:
            self.s = r_list
        
        try:
            # norm the unitvectors
            norms = N.sqrt((n_list**2).sum(axis=-1))
            assert N.all(norms > 0.)
        except AssertionError:
            raise ValueError("Unnormalizable orientation vector(s)")

        self.n = n_list / norms[:,N.newaxis]

        # assign other variables
        self.v = 0.
        self.P = P
        self.a = a
        self.l = l

        self.randState = randomUnitVector.RandomState(randSeed)
        self.noiseStDev = 0.

        if self.tumbles:
            self.t = 0.
        if self.sinks:
            self.F = F
        if self.forced:
            self.potentials = potentials.PotentialManager()
            
        return
    
    @staticmethod
    def _interp(lattice, r):
        return utils.interpArray(lattice, r)
    
    def move(self, lattice):
        """Updates the swimmers' positions using:

          Rdot = v(R) + Fn_/(6 pi eta a)

        where v(R) is the interpolated velocity at the position of the
        swimmer, a is the radius and n_ is the orientation.
        """
        rDot = N.zeros([self.num,3])   
        if self.swims:
            rDot += self.P * self.n 
            pass
        if self.sinks:
            rDot += self.F 
            pass
        if self.forced:
            rDot += self.G
            pass
        
        h = self.hydroRadius
        rDot *= (1./self.a - 1./h) / (6. * N.pi * self.eta)

        v = self._interp(lattice, self.r)
        if self.advects:
            rDot += v
            pass
        
        # rPlus = self.r
        nDot = N.zeros([self.num,3])   
        if self.rotates: 
            rMinus = self.r - self.n * self.l
            # from above, v = v(rPlus)
            nDot += v - self._interp(lattice, rMinus)
            # now nDot = v(rPlus) - v(rMinus)
            # so divide by l to get the rate
            nDot /= self.l
            pass

        self.applyMove(lattice, rDot)
        self.applyTurn(lattice, nDot)
        
        return
    
    def applyMove(self, lattice, rdots):
        """Moves the swimmer by an amount inc, respecting the periodic
        boundary conditions."""
        
        #lattice = self.L
        # set velocity
        self.v = rdots
        
         # new position
        r_ = self.r + rdots
       
       
        if self.walls:
            ind = N.where(r_[:,2] < 2)
            if len(ind[0]):
                z0 = self.r[ind][:, 2]
                #z1 = r_[ind][:, 2]
                # times to bottom is tau
                tau = (2. - z0) / rdots[ind][:, 2]
                # set intermediate posn
                r_[ind] = self.r[ind] + rdots[ind] * tau[:, N.newaxis]
                
                # randomize direction, to upper half sphere
                for i in ind[0]:
                    self.n[i] = self.randState.unitVector(3)
                    self.n[i, 2] = N.abs(self.n[i, 2])
                    continue
                
                # move v*(1 - tau) in the new direction
                r_[ind] += norm(rdots[ind])[:, N.newaxis] * (1. - tau[:, N.newaxis]) * self.n[ind]
                pass

            ind = N.where(r_[:,2] > lattice.nz -1)
            if len(ind[0]):
                z0 = self.r[ind][:, 2]
                #z1 = r_[ind][:, 2]
                # times to top is tau
                tau = (lattice.nz-1 - z0) / rdots[ind][:, 2]
                # set intermediate posn
                r_[ind] = self.r[ind] + rdots[ind] * tau[:, N.newaxis]
                
                # randomize direction, to lower half sphere
                for i in ind[0]:
                    self.n[i] = self.randState.unitVector(3)
                    self.n[i, 2] = -N.abs(self.n[i, 2])
                    continue
                
                # move v*(1 - tau) in the new direction
                r_[ind] += norm(rdots[ind])[:, N.newaxis] * (1. - tau[:, N.newaxis]) * self.n[ind]
                pass
            pass
        
        # Check that it's not moved off the edge of the lattice
        min = N.array([0.5, 0.5, 0.5])
        size = lattice.size
        max = size + 0.5
        
        for dim in range(3):
            ind = N.where(r_[:,dim] < min[dim])
            r_[ind, dim] += size[dim]
            ind = N.where(r_[:,dim] > max[dim])
            r_[ind, dim] -= size[dim]
            continue

        self.r = r_

        if self.tracks:
            # note that we DO NOT check for pbc
            # -- we want it to go over if needed
            self.s += rdots

            if self.walls:
                # make sure the unwrapped coords do respect the walls
                # by copying in the actual z
                self.s[2] = self.r[2]
        return
    
    def applyTurn(self, lattice, nDots):
        """Rotate the swimmers"""
        new = self.n + nDots
        norms = N.sqrt(N.sum(new**2, axis=-1))
        new /= norms[:, N.newaxis]
        self.n = new
        if self.tumbles:
            self.t = 0.
            rands = self.randState.rand(self.num)
            ind = N.where(rands < self.tumbleProb)
            for i in ind[0]:
                s = self[i]
                s.n = self.randState.unitVector(3)
                s.t = 1.
                continue
            
        return

    @staticmethod
    def _addForcesAtPositions(lattice, force, r):
        return utils.addForcesAtPositions(lattice, force, r)
        
    def addForce(self, lattice):
        # tail end
        # only force is the propulsion
        r = self.r - self.n * self.l
        force = -self.P * self.n
        self._addForcesAtPositions(lattice, force, r)

        # head end
        # have drag force equal to total force on particle
        force *= -1
        if self.sinks:
            force += self.F
            pass
        if self.forced:
            self.G[:] = 0.
            self.potentials.addForces(self.G, self.r, lattice=lattice)
            force += self.G
            pass
        
        self._addForcesAtPositions(lattice, force, self.r)
        
        return
    
    def respectsWalls(self):
        if self.walls:
            return 'z'
        else:
            return ''
    
    pass
