__registry = {}
def register(f):
    try:
        mfd = __registry[f.__module__]
    except KeyError:
        mfd = __registry[f.__module__] = {}
        pass
    
    mfd[f.func_name] = f
    return f
