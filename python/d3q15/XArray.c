#include "XArray.h"
#define NO_IMPORT_ARRAY
#include "DQarrayobject.h"

static PyObject *last = NULL;
static PyArrayObject *table = NULL;
static PyArrayObject *data = NULL;

static void XArray_updateCache(PyObject *self) {
  
  if (self == last) {
    /* we saw this one last time and since the table & number of cols etc.
       is fixed, we don't need to update the values stored. */
  } else {
    /* we haven't so extract some C values from the Python structures */
    data = (PyArrayObject *)PyObject_GetAttrString(self, "data");
    /* get the lookup table for columns */
    table = (PyArrayObject *)PyObject_GetAttrString(self, "cTable");
  }

}

int XArray_hasKey(PyObject *self, char k) {
  int nRows, nCols, i;
  char *array;
  
  XArray_updateCache(self);
  
  nRows = table->dimensions[0];
  nCols = table->dimensions[1];
  array = table->data;
  
  for (i=0; i<nRows; i++) {
    if (array[i*nCols + 0] == k) {
      /* found the right row */
      return 1;
    }
  }
  return 0;
}

char XArray_getStart(PyObject *self, char k) {
  int nRows, nCols, i;
  char *array;

  XArray_updateCache(self);
  
  nRows = table->dimensions[0];
  nCols = table->dimensions[1];
  array = table->data;
  
  for (i=0; i<nRows; i++) {
    if (array[i*nCols + 0] == k) {
      /* found the right row */
      return array[i*nCols + 1];
    }
  }
  PyErr_Format(PyExc_AttributeError, "XArray has no attribute '%c'", k);
  return XArray_ERROR;
}

char XArray_getLen(PyObject *self, char k) {
  int nRows, nCols, i;
  char *array;
  
  XArray_updateCache(self);
      
  nRows = table->dimensions[0];
  nCols = table->dimensions[1];
  array = table->data;
  
  for (i=0; i<nRows; i++) {
    if (array[i*nCols + 0] == k) {
      /* found the right row */
      return array[i*nCols + 1];
    }
  }
  PyErr_Format(PyExc_AttributeError, "XArray has no attribute '%c'", k);
  return XArray_ERROR;
}

double XArray_getValue(PyObject *self, int i, char key, int pos) {
  int nCols;
  char *array;
  
  XArray_updateCache(self);
  
  nCols = data->dimensions[1];
  array = data->data;
  
  return array[nCols*i + XArray_getStart(self, key) + pos];
}

void XArray_setValue(PyObject *self, int i, char key, int pos, double val) {
  int nCols;
  char *array;
  
  XArray_updateCache(self);
  
  nCols = data->dimensions[1];
  array = data->data;
  
  array[nCols*i + XArray_getStart(self, key) + pos] = val;
}
