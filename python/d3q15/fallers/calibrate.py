import d3q15
N = d3q15.N
import scipy
import cPickle

nIters = 100

def calibrateFromInstance(L, fallerClass):
    """Performs a calibration of the fitting parameter.
    L = the lattice to calibrate for
    fallerClass = the class of faller to calibrate for
    """
    latticeSize = (L.nx, L.ny, L.nz)
    return calibrate(L.__class__, latticeSize, L.tau_s/3.0, fallerClass, a, F_)

def calibrate(latticeClass, latticeSize, eta, fallerClass, posfile=None):
    """Performs a calibration of the fitting parameter.
    latticeClass = the class of lattice to calibrate for
    latticeSize = size of lattice
    eta = viscosity
    fallerClass = the class of faller to calibrate for
    """
    rng = N.random.RandomState()
    
    kT = 0.001
    
    (nx, ny, nz) = latticeSize
    tau_s = 3.0 * eta
    
    # no. steps to equilibriate the fluid
    equiliSteps = int((nx*ny*nz)**(2./3.) / eta)
    
    # gonna run for this many steps ONCE only and save the distributions
    L = d3q15.Lattice(nx,ny,nz, tau_s, tau_s)
    L.noise.temperature = kT
    L.noise.seed = rng.tomaxint()
    L.rho = N.ones(L.scalarFieldSize(), N.float)
    L.u = N.zeros(L.vectorFieldSize(), N.float)
    L.initBoundaryC('periodic')
    L.initForceC('none')
    L.initFromHydroVars()
    L.step(equiliSteps)
    
    fStore = L.f.copy()
    
    # Number of particles we can add if they are, on average, separated
    # by 4 lattice units (required for non-overlapping interpolation
    # regions
    nCells = nx*ny*nz
    nPart = max(nCells / 64, 1) # must be at least one!
    nRuns = nIters * nPart
    
    # How long to run for? Want to diffuse order 1 lattice unit to sample
    # over the grid. Use \Delta r ^2 ~ 2Dt = ____kT______ t
    #                                        3 \pi \eta a
    # With a = 1.5 (cos I know that's roughly right!)
    nSteps = int(3.*N.pi * eta * 1.5 / kT)
    
    posns = N.zeros((nSteps, nRuns, 3), N.float)
    
    for i in xrange(nIters):
        # construct initial particle positions
        initialPos = N.concatenate((rng.uniform(1,nx,size=(nPart,1)),
                                    rng.uniform(1,ny,size=(nPart,1)),
                                    rng.uniform(1,nz,size=(nPart,1))), axis=1)
        
        L = latticeClass(nx, ny, nz,tau_s, tau_s)
        L.initBoundaryC('periodic')
        L.add(fallerClass,
                     a=1., F=[0.,0.,0.], r_list=initialPos,
                     hydroRadius=1.1) 
        # turn off noise associated with subgridness of particle,
        L.things[0].noiseStDev = 0.0
        
        # copy in equilibriated fluid
        L.f[:] = fStore
        
        L.noise.temperature = kT
        # set a new random seed
        L.noise.seed = rng.tomaxint()
        
        L.updateHydroVars()
        
        # Going to integrate velocities ourselves to avoid PBC issues
        for j in xrange(nSteps):
            posns[j, nPart*i:nPart*(i+1), :] = \
                     posns[j-1, nPart*i:nPart*(i+1), :] + L.things[0].v
            L.step(1)
            continue
        
        continue
    
    if posfile is not None:
        cPickle.dump(posns, posfile, protocol=2)
        pass
    
    rSq = N.mean(N.sum(posns**2, axis=-1), axis=-1)
    t = N.arange(nSteps)
    coeffs = scipy.polyfit(t, rSq, 1)
    hydro = 2 * kT / (6 * N.pi * eta * coeffs[0])
    
    return hydro

