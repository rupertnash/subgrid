#ifndef XArray_H
#define XArray_H

#define XArray_ERROR 255

#include <Python.h>

int XArray_hasKey(PyObject *self, char k);
char XArray_getStart(PyObject *self, char k);
char XArray_getLen(PyObject *self, char k);
double XArray_getValue(PyObject *self, int i, char key, int pos);
void XArray_setValue(PyObject *self, int i, char key, int pos, double val);

#endif
