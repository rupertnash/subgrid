#ifndef Fallers_helpers_H
#define Fallers_helpers_H

#include <Python.h>
#include "d3q15.h"

/* static void FallerArray_addForce(PyObject *self, Lattice *lat); */
/* static PyObject *FallerArray_interp(Lattice *lat, PyObject *rArray); */
void FallerArray_move(PyObject *self, Lattice *lat);
void WalledFallerArray_move(PyObject *self, Lattice *lat);
void PDFallerArray_addPotentialDipoles(PyObject *self, Lattice *lat);

#endif
