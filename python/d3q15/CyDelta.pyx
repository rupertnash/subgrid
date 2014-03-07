"""Cython module implementing Peskin's regularization of Dirac's delta
function.

delta is a numpy Ufunc calculating the 1D regularization for every
element of its argument.

"""

cdef extern from "math.h":
    double fabs(double)
    double sqrt(double)

from numpy cimport npy_intp, NPY_DOUBLE, import_array
import_array()

    
cdef extern from "numpy/ufuncobject.h":
    void import_ufunc()
    ctypedef void (*PyUFuncGenericFunction) (char **, npy_intp *, npy_intp *, void *)
    object PyUFunc_FromFuncAndData(PyUFuncGenericFunction *, void **, char *, int, int, int, int, char *, char *, int)
    enum UFunc_macrosblah:
        PyUFunc_None
    

cdef double CyDelta_delta(double x):
    cdef double abs_x = fabs(x)
    cdef double root = -4. * x*x
    cdef double phi = -2.* abs_x
    
    if abs_x >= 2.0:
        return 0.
    
    if abs_x >= 1.0:
        root += 12. * abs_x - 7.
        phi += 5.
        phi -= sqrt(root)
    else:
        root += 4. * abs_x + 1
        phi += 3.
        phi += sqrt(root)
    
    return 0.125 * phi

cdef void delta_loop_1d(char **args, npy_intp *dims, npy_intp *steps, void *extra):
    cdef:
        npy_intp i
        npy_intp inpstep = steps[0]
        npy_intp outstep = steps[1]
        npy_intp n = dims[0]
        char *inp = args[0]
        char *out = args[1]
    
    for i in range(n):
        (<double *>out)[0] = CyDelta_delta((<double *>inp)[0])
        inp += inpstep
        out += outstep
    
    return

cdef PyUFuncGenericFunction delta_functions[1]
delta_functions[0] =  delta_loop_1d
cdef char delta_sigs[2]
delta_sigs[0] = NPY_DOUBLE
delta_sigs[1] = NPY_DOUBLE

import_ufunc()
delta = PyUFunc_FromFuncAndData(delta_functions,
                                NULL, delta_sigs,
                                1, 1, 1, PyUFunc_None, "delta",
                                "C implementation of Peskin's regularisation of the delta function.", 0)

def _delta(double x):
    return CyDelta_delta(x)
