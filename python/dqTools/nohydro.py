"""Some tools for creating Lattice objects which don't do the LB
simulation steps, and therefore don't have HI between objects.

Also contains shortcuts for SwimmerArrays in such situations.
"""
import d3q15.swimmers
import d3q15.wallswimmers
N = d3q15.N

class NoHydro(object):
    def step(self, n):
        for t in range(n):
            self.updateForce()
            self.time_step += 1

            if self.warmup:
                continue
            
	    for th in self.things:
		th.move(self)
		continue
	    
            continue
	return
    pass

class LatticeWithStuff(NoHydro, d3q15.LatticeWithStuff):
    pass

class SolidWalledLatticeWithStuff(NoHydro, d3q15.SolidWalledLatticeWithStuff):
    pass

class ExternalFlowSwimmerArray(d3q15.swimmers.SwimmerArray):
    @staticmethod
    def addForcesAtPositions(lattice, force, r):
        return
    
    pass

class ZeroFlowSwimmerArray(ExternalFlowSwimmerArray):
    @staticmethod
    def _interp(lattice, r):
        return N.zeros(r.shape)
    pass

SwimmerArray = ZeroFlowSwimmerArray

class ExternalFlowWallSwimmerArray(ExternalFlowSwimmerArray,
                                   d3q15.wallswimmers.WallSwimmerArray): pass

class ZeroFlowWallSwimmerArray(ZeroFlowSwimmerArray,
                               ExternalFlowWallSwimmerArray): pass
