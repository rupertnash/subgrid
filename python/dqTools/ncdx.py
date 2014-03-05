"""A simpler interface to netCDF files.

"""
from pycdf import *
import numpy as N
import os.path
CREATE = NC.CREATE|NC.WRITE|NC.BIT64_OFFSET

def get_or_def_dim(ncf, name, length):
    try:
        len = ncf.dimensions()[name]
        dim = ncf.inq_dimid(name)
    except KeyError:
        dim = ncf.def_dim(name, length)
    return dim

class Positions(object):
    pass

class Data(object):
    _default_rank=None
    def __init__(self, data, rank=None):
        if rank is None:
            rank = self._default_rank
            pass
        assert rank is not None

        assert data.ndim > rank
        
        self.rank = rank
        self.rankstr = ['', ', vector', ', matrix'][rank]
        
        itemshape = []
        shape = list(data.shape)
        
        for i in range(rank):
            itemshape.insert(0, shape.pop())
            continue
        
        self.shape = tuple(shape)
        self.itemshape = tuple(itemshape)
        self.size = N.multiply.reduce(shape)
        self.ndim = len(shape)
        
        self.data = data
        return

    def define(self, ncf, name):
        ncdims = []
        for i, n in enumerate(self.shape):
            ncdims.append(get_or_def_dim(ncf, name+'_values_n%d' % i, n))
            continue
        
        for i,n in enumerate(self.itemshape):
            ncdims.append(get_or_def_dim(ncf, name+'_values_r%d' % i, n))
            continue
        
        return ncf.def_var(name+'_values', NC.FLOAT, ncdims)

    def nameString(self, name):
        return name + self.rankstr
    
    def writeto(self, var):
        var.put(self.data.astype(N.float32))
        return
    
    pass

class ScalarData(Data):
    _default_rank = 0
    pass
class VectorData(Data):
    _default_rank = 1
    pass
class MatrixData(Data):
    _default_rank = 2
    pass

class IrregularPositions(Positions):
    def __init__(self, dim, positions):
        assert dim in [1,2,3]
        self.dim = dim
        
        positions = N.atleast_2d(positions)
        assert positions.ndim == 2
        assert positions.shape[1] == dim

        self.size = positions.shape[0]
        self.shape = (self.size,)
        self.ndim = 1
        
        self.data = self.positions = positions
        return
    
    def define(self, ncf, basename):
        sizeDim = get_or_def_dim(ncf, basename+'_values_n0', self.size)
        dimDim = get_or_def_dim(ncf, basename+'_dim', self.dim)
        
        return ncf.def_var(self.positionString(basename),
                           NC.FLOAT, (sizeDim, dimDim))

    @staticmethod
    def positionString(basename):
        return basename+'_locations'
    
    def writeto(self, var):
        var.put(self.positions.astype(N.float32))
        return
    
    pass

class RegularPositions(Positions):
    def __init__(self, axes):
        dim = len(axes)
        assert dim in [1,2,3]
        self.dim = dim
        self.axes = axes
        for i, ax in enumerate(axes):
            ax.dim = dim
            ax.aID = i
            continue
        
        return
    
    def define(self, ncf, basename):
        return [a.define(ncf, basename)
                for i, a in enumerate(self.axes)]
        
    def positionString(self, basename):
        ans = ''
        for i in range(self.dim):
            ans += self.axes[i].mkvarname(basename)
            
            if self.dim>1:
                ans += ', product'
                pass
            
            ans += self.axes[i].regular
            if i < (self.dim - 1):
                ans += '; '
            continue
        return ans
    
    def writeto(self, vars):
        for i in range(self.dim):
            self.axes[i].writeto(vars[i])
        return
    
    pass


class Axis(object):
    
    def mkvarname(self, basename):
        return basename+'_axis_%d' % self.aID
    
    def define(self, ncf, basename):
        dimDim = get_or_def_dim(ncf, basename+'_naxes', self.dim)
        sizeDim = get_or_def_dim(ncf, self.mkdimname(basename), self.n)
        return ncf.def_var(self.mkvarname(basename), NC.FLOAT, (sizeDim, dimDim))
    
    pass

class RegularAxis(Axis):
    regular = ', compact'
    def __init__(self, origin=0., delta=1., dim=None, aID=None):
        self.origin = origin
        self.delta = delta
        self.n = 2
        self.dim = dim
        self.aID = aID
        return
    
    def mkdimname(self, basename):
        return basename+'compact_dim'
    
    
    def writeto(self, var):
        pos = N.zeros((self.n, self.dim), dtype=N.float32)
        pos[0, self.aID] = self.origin
        pos[1, self.aID] = self.delta
        var.put(pos)
        return
    pass

class IrregularAxis(Axis):
    regular = ''
    def __init__(self, points, dim=None, aID=None):
        # has to be floats not doubles
        points = N.array(points, dtype=N.float32).squeeze()
        points = N.atleast_1d(points)
        
        assert points.ndim == 2
        assert len(points) >= 1
        
        self.n = points.shape[0]
        self.dim = points.shape[1]
        self.points = points
        self.aID = aID
        return
    
    def mkdimname(self, basename):
        return self.mkvarname(basename)+'_len'
    
    def writeto(self, var):
        var.put(self.points)
        return
    
    pass

class Field(object):
    
    def __init__(self, posns, **kwargs):
        self.name = kwargs.pop('name', 'data')
        self.positions = posns
        
        assert self.checkdict(kwargs, lambda key, val: isinstance(val, Data)),\
               "All keyword args must be Data istances"
        
        self.dataDict = kwargs
        return
    
    def writeto(self, ncf):
        #ncf = CDF(filename, mode=CREATE)
        
        ncf.definemode()
        pos = self.positions.define(ncf, self.name)

        datavars = {}
        for dname, d in self.dataDict.iteritems():
            fullvarname = "%s_%s" % (self.name, dname)
            dat = d.define(ncf, fullvarname)
            dat.field = d.nameString(fullvarname)
            dat.positions = self.positions.positionString(self.name)
            datavars[dname] = dat
            continue
        
        ncf.enddef()

        self.positions.writeto(pos)
        
        [d.writeto(datavars[dname]) for dname, d in self.dataDict.iteritems()]
        
        return
    
    @classmethod
    def checkdict(cls, d, condition):
        
        for key, val in d.iteritems():
            if condition(key, val):
                # we're ok here
                pass
            else:
                # test fails
                return False
            continue
        
        return True
    
    pass

class Unconnected(Field):
    def __init__(self, posns, **kwargs):
        assert isinstance(posns, IrregularPositions)
        
        Field.__init__(self, posns, **kwargs)
        
        assert self.checkdict(self.dataDict,
                              lambda key,val: val.ndim == 1), \
                              "data must be a list of positions"
        
        assert self.checkdict(self.dataDict,
                              lambda key,val: val.size == posns.size), \
                              "data & posns must have same size"
        
        return
    
    pass

class Connected(Field):
    def __init__(self, posns, **kwargs):
        assert isinstance(posns, RegularPositions)
        
        Field.__init__(self, posns, **kwargs)
        return


    pass

class LoadedData(object):
    def __init__(self, dname, vname, v):
        self.vname = vname
        self.dname = dname
        self.v = v
        
class LoadedField(object):
    def __init__(self, name):
        self.name = name
        self.positions = None
        self.dataDict = {}
        return
    
    def __getattribute__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            dd = object.__getattribute__(self,'dataDict')
            if name in dd:
                return dd[name]
            raise
        return
    
    pass

class LoadedFields(object):
    def __init__(self, ncfName):
        self.ncf = CDF(ncfName)
        vars = self.ncf.variables()
        dims = self.ncf.dimensions()
        self.fields = {}
        for vname, vinfo in vars.iteritems():
            v = self.ncf.var(vname)
            attr = v.attributes()
            if 'field' in attr:
                # it's a data field
                splitFieldAttr = attr['field'].split(',')
                try:
                    rankStr = splitFieldAttr[1].strip()
                except IndexError:
                    rankStr = 'scalar'
                    
                fname, dname  = splitFieldAttr[0].split('_')
                try:
                    f = self.fields[fname]
                except KeyError:
                    f = self.fields[fname] = LoadedField(fname)
                    pass
                
                f.dataDict[dname] = Data(v.get(),
                                         rank={'scalar': 0,
                                               'vector': 1,
                                               'matrix': 2}[rankStr])
                if f.positions is not None:
                    continue
                
                posdesc = attr['positions']
                if 'product' in posdesc:
                    # connected
                    axes = []
                    for dimdesc in posdesc.split(';'):
                        parts = dimdesc.split(',')
                        aID = parts[0].split('_')[2]
                        dimvar = self.ncf.var(parts[0].strip())
                        dimdat = dimvar.get()
                        if 'compact' in parts:
                            origin = dimdat[0,aID]
                            delta = dimdat[1,aID]
                            axes.append(RegularAxis(origin,delta))
                        else:
                            axes.append(IrregularAxis(dimdat))
                            pass
                        continue
                    pos = RegularPositions(axes)
                else:
                    # unconnected
                    posvar = self.ncf.var(posdesc.strip())
                    pos = IrregularPositions(dims[posvar.dimensions()[1]],
                                             posvar.get())
                    pass
                f.positions = pos
                pass
            continue
        return
    
    def __getattribute__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            f = object.__getattribute__(self,'fields')
            if name in f:
                return f[name]
            raise
        return

    pass

class Series(object):
    def __init__(self, dir, **kwargs):
        
        #mode='r', seriesvar='time', template='%.9d', clobber=False):
        try:
            func = {
            'r': self.openForRead,
            'w': self.create,
            'a': self.openForWrite
            }[kwargs.pop('mode', 'r')]
        except KeyError:
            raise ValueError('Invalid mode specified')
        
        func(dir, **kwargs)
        return
    
    def openForRead(self, dir, **kwargs):
        if not os.path.exists(dir):
            raise IOError('Series directory does not exist: %s' % dir)
        self.load(dir)
        return

    def openForWrite(self, dir, seriesvar='time', template='%.9d'):
        if not os.path.exists(dir):
            self.create(dir, seriesvar, template)
        else:
            self.load(dir)
        return
    
    def load(self, dir):
        self.dir = dir
        self.indexfile = os.path.join(dir, 'index')
        if not os.path.exists(self.indexfile):
            raise IOError('Series index does not exist: %s' % self.indexfile)
        
        header = file(self.indexfile).readline()
        if header[0] == '#':
            _, self.seriesvar = header[1:].split()
        else:
            self.seriesvar = 'time'
            pass
        self.index = [[el[0], eval(el[1])] for el in N.loadtxt(self.indexfile, object)]
        return
    
    def create(self, dir, seriesvar='time', template='%.9d'):
        if os.path.exists(dir):
            #clobber it
            import shutil
            shutil.rmtree(dir)
            pass
        
        self.dir = dir
        os.mkdir(dir)
        self.indexfile = os.path.join(dir, 'index')
        file(self.indexfile, 'w').write('# path %s\n' % seriesvar)
        self.seriesvar = seriesvar
        self.template = template
        self.index = []
        return
    
    def nextFile(self, var):
        filebase = (self.template % var) + '.ncdx'
        filename = os.path.join(self.dir, filebase)

        self.index.append([filebase, var])
        file(self.indexfile, 'a').write('%s %s\n' % (filebase, repr(var)))
        ncf = CDF(filename, mode=CREATE)
        return ncf
    
    def append(self, var, fields):
        ncf = self.nextFile(var)
        for f in fields:
            f.writeto(ncf)
            continue
        return

    def __len__(self):
        return len(self.index)
    
    def __getitem__(self, key):
        i = bisect_left(self.index, key)
        if self.index[i][1] == key:
            return LoadedFields(os.path.join(self.dir, self.index[i][0]))
        else:
            raise IndexError('Key "%s" not in index of "%s"' % (repr(key), self.dir))
        return
    
    def __iter__(self):
        for k in self.keys():
            yield self[k]
            
    def __contains__(self, key):
        i = bisect_left(self.index, key)
        if self.index[i][1] == key:
            return True
        else:
            return False
        return
    
    def keys(self):
        return [el[1] for el in self.index]
    pass

def bisect_left(a, x):
    lo = 0
    hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if a[mid][1] < x: lo = mid+1
        else: hi = mid
    return lo
def bisect_right(a, x):
    hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if x < a[mid][1]: hi = mid
        else: lo = mid+1
    return lo

# def scalarScatter(filename, values, positions, dim=3, name='data'):
#     Unconnected.write(filename, values, positions, dim=dim, rank=0, name=name)
# def vectorScatter(filename, values, positions, dim=3, name='data'):
#     Unconnected.write(filename, values, positions, dim=dim, rank=1, name=name)
# def matrixScatter(filename, values, positions, dim=3, name='data'):
#     Unconnected.write(filename, values, positions, dim=dim, rank=2, name=name)

