#ifndef utils_H
#define utils_H

#include <Python.h>
#include "d3q15.h"

void utils_addForceAtPosition_single(Lattice *lat,
				     double F[3],
				     double r[3]);

void utils_addForceAtPosition_array(Lattice *lat,
				    PyObject *FObj,
				    PyObject *rObj);

void utils_interp_single(Lattice *lat, double r[3], double v[3]);
void utils_interp_array(Lattice *lat, PyObject *rObj, PyObject * vObj);

#endif
