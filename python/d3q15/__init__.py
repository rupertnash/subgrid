"""This is a full-featured lattice Boltzmann implementation, using a
three-dimensional, fifteen-velocity model. The collision operator has
three relaxation times (shear, bulk and non-hydrodynamic, the last
being projected out at every timestep), as set out in Nash et
al. (DOI: 10.1103/PhysRevE.77.026709). Thermal fluctuations are also
incorporated with the method of Adhikari et al. (DOI:
10.1209/epl/i2004-10542-5).

Active and inactive particles are included through the swimmers and
fallers modules respectively. See my PhD thesis for full details.

d3q15.d3q15 is the SWIG proxy module
d3q15._d3q15 is the SWIG extension module

"""

from __future__ import absolute_import

import numpy as N
from . import d3q15, eigenvectors
import cPickle
import os.path

class Lattice(d3q15.Lattice):
    """Primary lattice Boltzmann class. Represents the computational
    domain, which must be a cuboid of fixed size. It is a further
    layer of wrapping around the SWIG proxy class which creates a more
    pythonic interface. 

    All quantities here are specified in lattice units \Delta x,
    \Delta t, etc..

    Indexing is zero-based but the fields all possess a halo of one
    lattice point's thickness, so the active parts are accessed with
    1..n inclusive.

    After construction, you MUST set the force update function with
    either:
     * setForceC -- uses a C-level function, must be one of those
                    specified in the documentation of setForceC
    or:
     * setForcePy -- any Python function object satistfying the 
                     criteria in the docstring of setForcePy

    Further you MUST set a boundary condition function (currently this
    must be in C, but could be extended similarly to the forcing),
    using initBoundaryC.
    
    """
    
    def __init__(self, *args):
        """Constructor. Arguments:
        int nx -- number of lattice points in the x-direction
        int ny -- number of lattice points in the y-direction
        int ny -- number of lattice points in the z-direction
        double tau_s -- relaxation time for the shear stress
        double tau_b -- relaxation time for the bulk stress
        
        The following all allow direct access to the LB fields,
        including the halos.
        f -- distributions
        rho -- density
        u -- velocity
        force -- force (actually acceleration density but should
        set rho==1)
        
        time_step -- the current time step
        
        The following are read-only computed properties:
        totalMass
        totalMomentum
        
        Note that the following attributes should be considered
        read-only but are in fact modifiable. They are provided here
        only for your reference, changing them will not affect the LB
        methods.
        
        size -- numpy array of the lattice size (nx, ny, nz)

        active -- an object that can be used to get the active
        (i.e. non-halo) parts of a Lattice's fields, like
        >>> lat.u[lat.active]

        eta -- the viscosity of the system
        
        w -- weights of each f_i
        
        xi -- velocity vectors
        
        Q -- kinetic projectors 

        complement -- complement of each velocity (xi[complement[i]]
                      == -xi[i])

        mm -- mode matrix that converts from the velocity basis to the
              moment basis (density, momentum, stress, ghosts, etc.;
              ordering as specified in libd3q15/eigh) (i.e. m_i = 
              mm_ij f_j)
        
        mmi -- inverse of mm
        
        norm -- normalizers for mm
        
        """
        d3q15.Lattice.__init__(self, *args)
        # now the basic object has been created, we can put a
        # pointer in the C struct to the PyObject
        self.__updateBackPointer__(self)
        
        # Create the lattice vectors, weights, kinetic projectors,
        # velocity<->mode matrix etc., as specified at compile time.
        eigenvectors.init(self)

        self.size = N.array((self.nx,
                             self.ny,
                             self.nz), dtype=int)
        
        self.active = (slice(1,self.nx+1),
                       slice(1,self.ny+1),
                       slice(1,self.nz+1))

        self.eta = self.tau_s/3.
        # some sensible default values for fields
        self.rho[:] = 1.0
        self.u[:] = 0.
        self.force[:] = 0.
        
        self.initFromHydroVars()
        return

    pass
    
class LatticeWithStuff(Lattice):
    """Implements a Lattice which contains things. See Lattice's
    docstring for full details.
    
    The things attribute is a list of the objects contained. This
    class will manage adding their contributions to the force to the
    fluid.

    You MUST NOT set the forcing function (this is taken care of) but
    you STILL MUST set the boundary condition.
    
    """
    def __init__(self, *args):
	Lattice.__init__(self, *args)
	self.things  = []
        # Note use of UNBOUND method so that the "self" argument isn't
        # automatically bound to the instance.
	self.setForcePy(self.__class__.updateForce)
	self.force[:] = 0.0
        self.warmup = False
	return
    
    def step(self):
        """Steps the system forward n time steps, including updating the
        positions of the Things.
        
        """
        Lattice.step(self)
        assert N.alltrue(N.isfinite(self.f_current.flat))
        self.updateHydroVars()
        assert N.alltrue(N.isfinite(self.u.flat))

        if self.warmup:
            return

        for th in self.things:
            th.move(self)
            continue
        return
    
    def add(self, thingClass, *args, **kwargs):
        """Calls the constructor/factory function supplied, and adds
        the resulting instance to the system. The instance must be a
        subclass of Thing and you must supply all necessary arguments
        after its class/factory.
        
        """
        thing = thingClass(self, *args, **kwargs)
        assert isinstance(thing, Thing)

        self.things.append(thing)
        thing.addForce(self)
	return thing
    
    def updateForce(self):
        """A method to be called from the C code to update the
        forcing onto the lattice.
        
        You may call this whenever it is neccessary to recompute the
        force field.
        """

        if self.warmup:
            return
        
	# zero the force first then add each contribution
	self.force[:] = 0.

        # The BodyForce, if present, must be done last
	bf = None
        
	for th in self.things:
            if isinstance(th, BodyForce):
                # It's the BodyForce, so store for later
                bf = th
            else:
                # Add the force contribution normally
                th.addForce(self)
	    continue

        if bf is not None:
            # Add the BodyForce now
            bf.addForce(self)
            
	return
    

    def __getstate__(self):
        """Enables pickling - creates the actual data to be
        serialized."""
	picdict = Lattice.__getstate__(self)
	# Remove rho, u and force field since they can be recalculated
	# from other data.
	del picdict['fields']['rho']
	del picdict['fields']['u']
	del picdict['fields']['force']
	picdict['things'] = self.things
        picdict['warmup'] = self.warmup
	return picdict
    
    def __setstate__(self, picdict):
        """Restores from the pickled data."""
	try:
	    # Since we've removed rho, u & force from the 'fields' 
	    # sub-dict, KeyError will occur. Ignore it
	    Lattice.__setstate__(self, picdict)
	except KeyError:
	    pass
	
	self.things = picdict['things']
        # initially False to ensure the force is updated onece
        self.warmup = False
	# But now need to reconstruct the fields
	self.updateForce()
	self.updateHydroVars()
        # Now set warmup flag to correct state.
        try:
            self.warmup = picdict['warmup']
        except KeyError:
            self.warmup = False
                
	return
    
    pass


class Thing(object):
    """A 'virtual' class for things to be placed in the Lattice. Your
    classes must override the addForce and move methods.

    """
    def __init__(self, *args, **kwargs):
        return
    
    def addForce(self, lattice):
	raise NotImplementedError("Subclass must override this method!")
    
    def move(self, lattice):
	raise NotImplementedError("Subclass must override this method!")
    
    pass

class SolidWalledLatticeWithStuff(LatticeWithStuff):
    """Implements a Lattice which contains things and has stationary,
    solid walls in the z-direction. See LatticeWithStuff's & Lattice's
    docstrings for full details.
    
    Added Things must have a method respectsWalls that returns 'z'.
    
    You MUST NOT set the forcing function (this is taken care of).
    You MUST NOT set the boundary condition function (this is taken care of).
    
    """
    def __init__(self, *args):
	LatticeWithStuff.__init__(self, *args)
	self.initBoundaryC('noslip')
	return
    
    def add(self, thingClass, *args, **kwargs):
	th = LatticeWithStuff.add(self, thingClass, *args, **kwargs)
        try:
            assert th.respectsWalls() == 'z'
        except AssertionError:
            del self.things[-1]
            raise
        
	return th
    
    pass

	
class BodyForce(Thing):
    """A special Thing that ensures there is no net force on the
    LatticeWithStuff. It's special-cased by the LatticeWithStuff
    class. You can use this as a superclass of anything that requires
    this behaviour, although there can be only one and this is NOT
    CHECKED!
    
    """
    
    def addForce(self, lat):
        totalForce = lat.force[lat.active].sum(axis=0).sum(axis=0).sum(axis=0)
        nSites = lat.nx * lat.ny * lat.nz
        
        forcePerSite = -totalForce / nSites
        
        lat.force += N.reshape(forcePerSite, (1,1,1,3))
        return

    def move(self, lattice):
        return
    
    pass

    
class XArray(Thing):
    """This is a base class for fixed-size arrays of simple objects
    that can be described by a series of one-dimensional vectors. Each
    vector attribute is keyed by a single character. There is a
    corresponding C-API for accessing these properties (see XArray.h).
    
    The attributes must be created at initialization time by calls to
    tableAppend, typically in the subclass's __init__. The method MUST
    set the "num" attribute to the number of objects to be stored in
    the array, and MUST then call the superclass.
    
    After initialization, the attributes can be simply accessed by the
    usual attribute reference (object.attributename) and single items
    from the array retrieved by indexing (object[i]).

    Example
    -------
    from d3q15 import XArray
    import numpy as N
    
    class ChargeArray(XArray): 
        def __init__(self, L, points=None, charges=None):
            assert isinstance(points, N.ndarray)
            assert len(points.shape) == 2
            assert points.shape[1] ==3
            
            assert isinstance(charges, N.ndarray)
            assert len(charges.shape) == 1
            assert points.shape[0] == charges.shape[0]
            
            self.num = points.shape[0]
            # position
            self.tableAppend('r', 3)
            # charge
            self.tableAppend('q', 1)
            XArray.__init__(self, L)
            
            self.r = points
            self.q = charges
            
    lat = initLB()
    initialPos = N.array([[3., 4., 5.],
                          [1.5, 3., 10.],
                          [2.,3.,4.]])
    charge = N.ones(3)
    
    pa = PointArray(lat, points=initialPos)
    
    # pa.q is a numpy array of all charges
    print pa.q
    # [ 1.  1.  1.]
    point = pa[1]
    print point.r
    # [ 1.5  3.  10.]
    
    pa.q = -1
    print point.q
    # [ -1.]
    
    """
    
    def __init__(self, L, **kwargs):
        """Initializes the data array and the lookup tables
        required."""
        if not hasattr(self, 'table'):
            self.table = {}
            pass
        #self.L = L
        self._initCTable()
        self._initArrays()
        return
    
    def tableAppend(self, key, length):
        """Add a new attribute to the XArray.
        
        key -- a single character which is name of the attribute
        length -- an integer specifying the length of the 1D vector
        """
        assert isinstance(key, str) and len(key) == 1, \
               "All keys must be length 1 strings"
        assert key not in self.table, "No duplicate keys"
        
        curLen = 0
        for k in self.table.keys():
            curLen += self.table[k][1]
            continue
        
        self.table[key] = (curLen, length)
        return

    def _initCTable(self):
        """Creates a simple table that is accessed with the C API."""
        t = self.table
        nColumns = 0
        for k in t.keys():
            assert isinstance(k, str) and len(k) == 1, \
                   "All keys must be length 1 strings"
            nColumns += t[k][1]
            continue
        self.nColumns = nColumns

        # For the c code, we limit the number of columns to 256 by
        # using char data. Pad with an extra zero to ensure each
        # row is 32 bits wide.
        self.cTable = N.zeros((len(t), 4),
                             dtype=N.byte)
        
        # The keys are one character strings
        keys = t.keys()
        for i in range(len(t)):
            k = keys[i]
            self.cTable[i, :] = [N.fromstring(k, dtype=N.byte)[0],
                                 t[k][0],
                                t[k][1],
                                0]
            continue
        return

    def _get_nCols(self):
        table = self.table
        nColumns = 0
        for k in table.keys():
            assert isinstance(k, str) and len(k) == 1, \
                   "All keys must be length 1 strings"
            nColumns += table[k][1]
            continue
        return nColumns
    
    def _initArrays(self):
        self.data = N.zeros((self.num, self._get_nCols()), N.float)
        return

    def __getitem__(self, i):
        return self.singleClass(self, i)
    
    def __iter__(self):
        for i in range(self.num):
            yield self[i]

    def __getattr__(self, name):
        try:
            lwr,ln = object.__getattribute__(self, 'table')[name]
            return object.__getattribute__(self, 'data')[:, lwr:lwr+ln]
        except AttributeError, KeyError:
            raise AttributeError("'%s' object has no attribute '%s'" % (self.__class__.__name__, name))
        
    def __setattr__(self, name, value):
        try:
            table = object.__getattribute__(self, 'table')
            if name in table.keys():
                lwr,ln = table[name]
                object.__getattribute__(self, 'data')[:, lwr:lwr+ln] = value
                return
            
        except AttributeError:
            pass
        
        object.__setattr__(self, name, value)

    pass

class X(object):
    """A single entry from an XArray. These objects can only exist as
    a effectively a row from the table that is an XArray instance.

    You should probably never directly use this class; see XArray
    docstring for an example of use.
    
    """
    def __init__(self, container, i):
        object.__setattr__(self, 'container', container)
        object.__setattr__(self, '_data', container.data[i])
        return
    
    def __getattr__(self, name):
        try:
            lwr,ln = object.__getattribute__(self, 'container').table[name]
            return object.__getattribute__(self, '_data')[lwr:lwr+ln]
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'" % (self.__class__.__name__, name))
        
    def __setattr__(self, name, value):
        # must use self.__dict__ otherwise everything can get very confused
        # since __setattr__ is always called
        # Note hasattr implicitly calls __getattr__
        
        try:
            lwr,ln = self.container.table[name]
            self._data[lwr:lwr+ln] = value
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'" % (self.__class__.__name__, name))
        

    pass

XArray.singleClass = X
