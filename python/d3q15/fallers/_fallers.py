import d3q15
from d3q15.CDelta import delta
import d3q15.utils
import _helpers # d3q15.fallers._helpers as _helpers

N = d3q15.N

        
class FallerArray(d3q15.XArray):

    def __init__(self, L, **kwargs):
        """Create an array of Faller objects.

        Arguments:
          L: the Lattice object on which the Fallers are falling

        Keyword arguments:
          a: the radius of the particles
          r_list: the initial positions of the particles
          F: the force vector that is acting on the particles
          hydroRadius: the lattice hydrodynamic radius, use the
          module fallers.calibrate.calibrate to determine its value

        Optional keyword arguments:
          randSeed: seed for PRNG, seeded from /dev/urandom if absent
          walls: flag indicating whether the Lattice has walls
          pd: flag indicating whether to use the potential dipole term
          static: flag indicating whether to keep the particles fixed
        """
        try:
            a = kwargs['a']
            F = kwargs['F']
            r_list = kwargs['r_list']

        except KeyError, err:
            raise TypeError("Must specify %s" % err)

        try:
            hydroRadius = kwargs['hydroRadius']
        except KeyError:
            raise TypeError("Must specify hydroRadius. Use fallers.calibrate.calibrate to determine its value")

        try:
            randSeed = kwargs['randSeed']
        except KeyError:
            randSeed = 0

        for flag in ['walls', 'pd', 'static']:
            try:
                if kwargs[flag]:
                    self.__setattr__(flag, True)
                else:
                    self.__setattr__(flag, False)
            except KeyError:
                self.__setattr__(flag, False)
            continue
        
        self.eta = L.tau_s / 3.
        self.num = len(r_list)

        # set up fields in data array
        # we only create a table if one doesn't exist already
        # (we might be called by a subclass that created some already)
        if not hasattr(self, 'table'):
            self.table = {}
        self.tableAppend('r',3)
        self.tableAppend('v',3)
        self.tableAppend('F',3)
        self.tableAppend('a',1)
        d3q15.XArray.__init__(self, L)
            
        if not (N.alltrue(r_list[:,0] > 0.5) and
                N.alltrue(r_list[:,0] <= L.nx+0.5)):
            raise ValueError("xcoordinate out of range")
        if not (N.alltrue(r_list[:,1] > 0.5) and
                N.alltrue(r_list[:,1] <= L.ny+0.5)):
            raise ValueError("ycoordinate out of range")
        if not (N.alltrue(r_list[:,2] > 0.5) and
                N.alltrue(r_list[:,2] <= L.ny+0.5)):
            raise ValueError("zcoordinate out of range")

        if self.walls:
            if not (N.alltrue(r_list[:,2] > 2) and
                    N.alltrue(r_list[:,2] <= L.nz - 1)):
                raise ValueError("z-coordinate out of range due to walls")
        
        self.r = r_list
        self.v = 0.
        self.F = F
        self.a = a
        
        assert hydroRadius > a
        self.hydroRadius = hydroRadius
        self.rng = N.random.RandomState(randSeed)
        self.noiseStDev = N.sqrt(
            L.noise.temperature/(9.*N.pi*self.eta) * (1./a - 1./hydroRadius)
            )
        
        return

    def move(self, lattice):
        if not self.static:
            if self.walls:
                _helpers.WalledFallerArray_move(self, lattice)
            else:
                _helpers.FallerArray_move(self, lattice)
        
        return

    def addForce(self, lattice):
        d3q15.utils.addForcesAtPositions(lattice, self.F, self.r)
        
        if self.pd:
            _helpers.PDFallerArray_addPotentialDipoles(self, lattice)
        return
    
    def noise(self, n):
        """Work out the noise required to make the particle diffuse correctly
        D = kT/(6 \pi \eta h)
        """
        return self.noiseStDev * self.rng.normal(size=3)
    
    def interpOne(self, lattice, i):
        """Delta function interpolation:
        v(R) = \int dr v(r) \delta (r - R)
        where we approximate the integral with a sum over the support of
        Peskin's delta function.
        """
        R = self[i].r
        v = d3q15.utils.interpOne(lattice, R)
        #assert N.alltrue(N.isfinite(v))
        return v
    
    def interpAll(self, lattice):
        return d3q15.utils.interpArray(lattice, self.r)
    
    def moveOne(self, lattice, i, inc):
        """Moves the faller by an amount inc, respecting the periodic
        boundary conditions."""
        # set the velocity
        self[i].v = inc
        size = N.array((lattice.nx, lattice.ny, lattice.nz))
        
        # Calculate the new position
        r_ = self[i].r + inc
        
        if self.walls and r_[2] < 2.:
            # it has passed thru the bottom allowed position, 
            # so set its force to zero, so it won't affect things anymore
            self.F = 0.
            pass
        
        # Now check that it's not moved off the edge of the lattice
        # and wrap it round if it has.
        ind = N.where(r_<0.5)
        r_[ind] += size[ind]
        ind = N.where(r_>size+0.5)
        r_[ind] -= size[ind]
        
        # store the new position
        self[i].r = r_
        return
    
    def respectsWalls(self):
        if self.walls:
            return 'z'
        else:
            return ''
        
    pass

    
def test():
    size = 10
    tau = 0.5
    eta = tau/3.

    posns = N.array([[1.,1.,1.],
                     [10,10,10]])
    a = 0.05
    hydroRadius = 1.5
    F = N.array([0,0,-1e-4])

    v_s = -F[-1]/(6*N.pi*eta*a)

    print "Basic faller test"
    lat = d3q15.LatticeWithStuff(size, size, size, tau,tau)
    lat.initBoundaryC('periodic')
    lat.add(d3q15.BodyForce)
    lat.add(FallerArray, r_list=posns, a=a, F=F, hydroRadius=hydroRadius)
    lat.initFromHydroVars()
    try:
        lat.step(int(1/v_s))
    except KeyboardInterrupt:
        print "Skipping..."
        pass

    print "PD faller test"
    lat = d3q15.LatticeWithStuff(size, size, size, tau,tau)
    lat.initBoundaryC('periodic')
    lat.add(d3q15.BodyForce)
    lat.add(FallerArray, r_list=posns, a=a, F=F, hydroRadius=hydroRadius, pd=True)
    lat.initFromHydroVars()
    try:
        lat.step(int(1/v_s))
        print "Heights should be approx: " +str(posns[:,2] - 1)
        print "Are: " + str(lat.things[-1].r[:,2])
    except KeyboardInterrupt:
        print "Skipping..."
        pass

    print "Static faller test"
    lat = d3q15.LatticeWithStuff(size, size, size, tau,tau)
    lat.initBoundaryC('periodic')
    lat.add(d3q15.BodyForce)
    lat.add(FallerArray, r_list=posns, a=a, F=F, hydroRadius=hydroRadius, static=True)
    lat.initFromHydroVars()
    try:
        lat.step(int(1/v_s))
    except KeyboardInterrupt:
        print "Skipping..."
        pass

    print "PD and static faller test"
    lat = d3q15.LatticeWithStuff(size, size, size, tau,tau)
    lat.initBoundaryC('periodic')
    lat.add(d3q15.BodyForce)
    lat.add(FallerArray, r_list=posns, a=a, F=F, hydroRadius=hydroRadius, pd=True, static=True)
    lat.initFromHydroVars()
    try:
        lat.step(int(1/v_s))
    except KeyboardInterrupt:
        print "Skipping..."
        pass

    posns = N.array([[1.,1.,2.5],
                     [10,10,9.]])
    print "Walled faller test"
    lat = d3q15.SolidWalledLatticeWithStuff(size, size, size, tau,tau)
    lat.add(FallerArray, r_list=posns, a=a, F=F, hydroRadius=hydroRadius, walls=True)
    lat.initFromHydroVars()
    try:
        lat.step(int(1/v_s))
    except KeyboardInterrupt:
        print "Skipping..."
        pass

