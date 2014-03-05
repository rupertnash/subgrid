"""Simple, testing module to add a list of delta-function forces to a lattice.


"""
from __future__ import absolute_import

from math import ceil
from .delta import delta
import numpy as N

__deltaList = None

def forceInit(lat, deltas):
    """Initialises the forcing.

    Note that this stores a reference to the list of delta functions, so
    any changes to that list will be reflected here.

    deltas must be a container type. Each element must be a pair of
    three-vectors: (position, force) = ((x,y,z), (Fx, Fy, Fz)).
    
    """
    global __deltaList
    __deltaList = deltas
    forceUpdate(lat)


def forceUpdate(lat):
    """Updates the force on the lattice, using the list of delta
    functions, overwrites any existing forces applied.
    """
    global __deltaList

    # zero the force field
    force = N.zeros((lat.nx+2, lat.ny+2, lat.nz+2, 3), N.float)
    for d in __deltaList:
        x = d[0]
        F = d[1]
        
        # Check x is in range
        if not(x[0] >= 0.0 and x[0] < lat.nx):
            raise ValueError("x coordinate out of range for delta " + str(d))
        if not(x[1] >= 0.0 and x[1] < lat.ny):
            raise ValueError("y coordinate out of range for delta " + str(d))
        if not(x[2] >= 0.0 and x[2] < lat.nz):
            raise ValueError("z coordinate out of range for delta " + str(d))

        xpoints = range(int(ceil(x[0]-2)), int(ceil(x[0]+2)))
        xrange = [(i-1)%lat.nx+1 for i in xpoints]
        ypoints = range(int(ceil(x[1]-2)), int(ceil(x[1]+2)))
        yrange = [(i-1)%lat.ny+1 for i in ypoints]
        zpoints = range(int(ceil(x[2]-2)), int(ceil(x[2]+2)))
        zrange = [(i-1)%lat.nz+1 for i in zpoints]

        for i in range(len(xrange)):
            for j in range(len(yrange)):
                for k in range(len(zrange)):
                    delta3D = delta(x[0]-xpoints[i]) * delta(x[1]-ypoints[j]) * delta(x[2]-zpoints[k])
                    for dim in range(3):
                        force[xrange[i],yrange[j],zrange[k], dim] += delta3D * F[dim]

        lat.force = force

        
