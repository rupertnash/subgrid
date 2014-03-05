"""A simple library for serializing data to disk and allowing random
access.

Probably best not to use this-- either straightforward pickling or a
"proper" scientific format like netCDF or HDF.

I don't think any subgrid runs use this, but I can't be sure.

"""
import cPickle
import os.path

# i.e. 0.5GB
MAXFILESIZE = 2**29
VERSION = 3

def new(name):
    """Factory function that creates a new RecordDB with filename name."""
    rdb = RecordDB()
    rdb.new(name)
    return rdb

def open(name, readonly=False):
    """Factory function that opens an existing RecordDB with filename name."""
    ver = _versionCheck(name)
    assert ver <= VERSION
    if ver == 1:
        rdb = RecordDB_V1()
    elif ver == 2:
        rdb = RecordDB_V2()
    else:
        rdb = RecordDB()
        
    rdb.open(name, readonly)
    return rdb

def _versionCheck(name):
    name = os.path.abspath(name)
    assert os.path.exists(name)
    assert os.path.isdir(name)
    
    v = 1
    for f in os.listdir(name):
        (base,ext) = os.path.splitext(f)
        if base == 'version':
            v = int(ext[1:])
            break
        continue
    
    return v


class RecordDB_V3(object):
    """Class that wraps an indexed store of records.
    Slightly altered to append only, making it suitable for use with
    auto-checkpointing!"""
    
    def __init__(self):
        self.version = 3
        self.readonly = False
        return
    
    def new(self, name):
        """Initialises a new RDB."""
        name = os.path.abspath(name)
        assert not os.path.exists(name)
        self.name = name
        os.mkdir(name)
        file(os.path.join(name, 'version.%d' % self.version), 'w').close()
                
        self.map = []
        self.mapFN = os.path.join(name, 'map.dat')
        file(self.mapFN, 'w').close()
		
        return

    def open(self, name, readonly):
        """Opens an existing RDB."""
        name = os.path.abspath(name)
        assert os.path.exists(name)
        assert os.path.isdir(name)
        self.name = name
        self.mapFN = os.path.join(self.name, 'map.dat')
        assert os.path.exists(self.mapFN)
        assert os.path.isfile(self.mapFN)
        
        self.loadMap()
        self.readonly = readonly
        
        return

    def __len__(self):
        return len(self.map)

    def __getitem__(self, key):
        if type(key) is int:
            (fn, posn, tag) = self.map[key]
        else:
            found = []
            for entry in self.map:
                (fn, posn, tag) = entry
                if tag == key:
                    found.append(entry)
                    pass
                continue
            if len(found) == 0:
                raise KeyError(key)
            if len(found) > 1:
                raise KeyError("Duplicate keys cannot be searched for!")
            (fn, posn, tag) = found[0]
            
                
        f = file(os.path.join(self.name,fn), 'rb')
        try:
            f.seek(posn)
            p = cPickle.Unpickler(f)
            return p.load()
        finally:
            f.close()
    
    def append(self, item, tag=None):
        """Appends item to the DB."""
        assert self.readonly == False
        
        fn = self.chooseFile(item)
        f = file(os.path.join(self.name,fn), 'ab')
        try:
            posn = f.tell()
            p = cPickle.Pickler(f, protocol=2)
            p.dump(item)
        finally:
            f.close()
        # update the map
        self.map.append((fn, posn, tag))
        
        mapfile = file(self.mapFN, 'ab')
        cPickle.dump((fn, posn, tag), mapfile, protocol=2)
        mapfile.close()
        return

    def __setitem__(self, key, value):
        raise NotImplementedError
    
    def loadMap(self):
        """Loads the map from the mapfile."""
        f = file(self.mapFN, 'rb')
        self.map = []
        while True:
            try:
                self.map.append(cPickle.load(f))
            except EOFError:
                break
        
        f.close()
        return
    
    def chooseFile(self, item):
        """Selects the filename to which to write the next object."""
        # get the last file
        if len(self.map) == 0:
            fn = "000000.dat"
        else:
            (fn, posn, tag) = self.map[-1]
            sizeB = os.path.getsize(os.path.join(self.name,fn))
            if sizeB > MAXFILESIZE:
                ind = int(fn[:6]) + 1
                assert ind <= 999999
                fn = "%.6d.dat" % ind
            
        return fn
    
class RecordDB_V2(object):
    """Class that wraps an indexed store of records."""

    def __init__(self):
        self.readonly = False
        return
    
    def new(self, name):
        """Initialises a new RDB."""
        name = os.path.abspath(name)
        assert not os.path.exists(name)
        self.name = name
        os.mkdir(name)
        file(os.path.join(name, 'version.%d' % 2), 'w').close()
        self.version = 2
        
        self.map = []
        self.mapFN = os.path.join(name, 'map.dat')
        self.dumpMap()
        
        return

    def open(self, name, readonly):
        """Opens an existing RDB."""
        name = os.path.abspath(name)
        assert os.path.exists(name)
        assert os.path.isdir(name)
        self.name = name
        self.mapFN = os.path.join(self.name, 'map.dat')
        assert os.path.exists(self.mapFN)
        assert os.path.isfile(self.mapFN)
        
        self.loadMap()
        self.readonly = readonly
        
        return

    def __len__(self):
        return len(self.map)

    def __getitem__(self, key):
        if type(key) is int:
            (fn, posn, tag) = self.map[key]
        else:
            found = []
            for entry in self.map:
                (fn, posn, tag) = entry
                if tag == key:
                    found.append(entry)
                    pass
                continue
            if len(found) == 0:
                raise KeyError(key)
            if len(found) > 1:
                raise KeyError("Duplicate keys cannot be searched for!")
            (fn, posn, tag) = found[0]
            
                
        f = file(os.path.join(self.name,fn), 'rb')
        try:
            f.seek(posn)
            p = cPickle.Unpickler(f)
            return p.load()
        finally:
            f.close()
    
    def append(self, item, tag=None):
        """Appends item to the DB."""
        assert self.readonly == False
        
        fn = self.chooseFile(item)
        f = file(os.path.join(self.name,fn), 'ab')
        try:
            posn = f.tell()
            p = cPickle.Pickler(f, protocol=2)
            p.dump(item)
        finally:
            f.close()
        # update the map
        self.map.append((fn, posn, tag))
        self.dumpMap()
        return

    def __setitem__(self, key, value):
        raise NotImplementedError

    def dumpMap(self):
        """Writes the map to the mapfile."""
        f = file(self.mapFN, 'wb')
        try:
            p = cPickle.Pickler(f, protocol=2)
            p.dump(self.map)
        finally:
            f.close()
        return
    
    def loadMap(self):
        """Loads the map from the mapfile."""
        f = file(self.mapFN, 'rb')
        try:
            p = cPickle.Unpickler(f)
            self.map = p.load()
        finally:
            f.close()
        return
    
    def chooseFile(self, item):
        """Selects the filename to which to write the next object."""
        # get the last file
        if len(self.map) == 0:
            fn = "000000.dat"
        else:
            (fn, posn, tag) = self.map[-1]
            sizeB = os.path.getsize(os.path.join(self.name,fn))
            if sizeB > MAXFILESIZE:
                ind = int(fn[:6]) + 1
                assert ind <= 999999
                fn = "%.6d.dat" % ind
            
        return fn
        
            
            
class RecordDB_V1(object):
    """Class that wraps an indexed store of records."""

    def __init__(self):
        self.readonly = False
        return
    
    def new(self, name):
        """Initialises a new RDB."""
        name = os.path.abspath(name)
        assert not os.path.exists(name)
        self.name = name
        os.mkdir(name)
        
        self.map = []
        self.mapFN = os.path.join(name, 'map.dat')
        self.dumpMap()
        
        return

    def open(self, name, readonly):
        """Opens an existing RDB."""
        name = os.path.abspath(name)
        assert os.path.exists(name)
        assert os.path.isdir(name)
        self.name = name
        self.mapFN = os.path.join(self.name, 'map.dat')
        assert os.path.exists(self.mapFN)
        assert os.path.isfile(self.mapFN)
        
        self.loadMap()
        self.readonly = readonly
        
        return

    def __len__(self):
        return len(self.map)

    def __getitem__(self, key):

        (fn, posn) = self.map[key]
        f = file(os.path.join(self.name,fn), 'rb')
        try:
            f.seek(posn)
            p = cPickle.Unpickler(f)
            return p.load()
        finally:
            f.close()
    
    def append(self, item):
        """Appends item to the DB."""
        assert self.readonly == False
        
        fn = self.chooseFile(item)
        f = file(os.path.join(self.name,fn), 'ab')
        try:
            posn = f.tell()
            p = cPickle.Pickler(f, protocol=2)
            p.dump(item)
        finally:
            f.close()
        # update the map
        self.map.append((fn, posn))
        self.dumpMap()
        return

    def __setitem__(self, key, value):
        raise NotImplementedError

    def dumpMap(self):
        """Writes the map to the mapfile."""
        f = file(self.mapFN, 'wb')
        try:
            p = cPickle.Pickler(f, protocol=2)
            p.dump(self.map)
        finally:
            f.close()
        return
    
    def loadMap(self):
        """Loads the map from the mapfile."""
        f = file(self.mapFN, 'rb')
        try:
            p = cPickle.Unpickler(f)
            self.map = p.load()
        finally:
            f.close()
        return
    
    def chooseFile(self, item):
        """Selects the filename to which to write the next object."""
        # get the last file
        if len(self.map) == 0:
            fn = "000000.dat"
        else:
            (fn, posn) = self.map[-1]
            sizeB = os.path.getsize(os.path.join(self.name,fn))
            if sizeB > MAXFILESIZE:
                ind = int(fn[:6]) + 1
                assert ind <= 999999
                fn = "%.6d.dat" % ind
            
        return fn
    
    pass

RecordDB = RecordDB_V3
