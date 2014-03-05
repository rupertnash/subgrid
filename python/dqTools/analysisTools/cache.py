import os.path
import cPickle
from functools import wraps

def cache(filename, *argPatterns, **kwargPatterns):
    """Given a name of a file, returns a decorator which will cache
    the results of a function call (of a function with a first arg
    of "base") to "base/filename"
    """
    def decorator(f):
        # this is the actual decorator
        # f is the function to be modified
        
        @wraps(f)
        def cacher(base, *args, **kwargs):
            # create the full path to the cache file
            cachefile = filename
            if len(argPatterns):
                cachefile += '.'+'.'.join([a%args[i] 
                                           for i,a in enumerate(argPatterns)])
                pass
            
            keys = kwargPatterns.keys()
            if len(keys):
                keys.sort()
                cachefile += '.'+'.'.join([kwargPatterns[k]%kwargs[k] 
                                           for k in keys])
                pass
            
            cacheFN = os.path.join(base, cachefile)
            
            if os.path.exists(cacheFN):
                ans = cPickle.load(file(cacheFN))
            else:
                ans = f(base, *args, **kwargs)
                cPickle.dump(ans,
                             file(cacheFN, 'wb'),
                             protocol=2)
                pass
            return ans
        
        # return the wrapped function
        return cacher

    return decorator
