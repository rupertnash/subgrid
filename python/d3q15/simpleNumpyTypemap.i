/* -*- mode: c; -*- */
%header %{
#include "DQarrayobject.h"
%}

%fragment("getOrMakeArrayOfDoubles", "wrapper") {
  /* helper function */
  static double *getOrMakeArrayOfDoubles(PyObject *input, int size, double *temp) {
    int i;
    PyObject *o = NULL;
    
    if (!PyArray_Check(input)) {
      if (!PySequence_Check(input)) {
	PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
	return NULL;
      }
      
      if (PyObject_Length(input) != size) {
	PyErr_Format(PyExc_ValueError,
		     "Expecting a sequence with %d elements", size);
	return NULL;
      }
      
      for (i=0; i<size; i++) {
	o = PySequence_GetItem(input, i);
	if (!PyFloat_Check(o)) {
	  Py_XDECREF(o);
	  PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats");
	  return NULL;
	}
	temp[i] = PyFloat_AsDouble(o);
	Py_DECREF(o);
      }
      return temp;
      
    } else {
      PyArrayObject *arrInput = (PyArrayObject *)input;
      int nelem = 1;
      /* check array of doubles */
      if (PyArray_DTYPE(arrInput)->type_num != NPY_DOUBLE) {
	PyErr_SetString(PyExc_TypeError, "Expecting an array of doubles");
	return NULL;
      }
      
      npy_intp* dims = PyArray_DIMS(arrInput);
      for (i=0; i<PyArray_NDIM(arrInput); i++) {
	nelem *= dims[i];
      }
      if (nelem!= size) {
	PyErr_Format(PyExc_ValueError,
		     "Expecting a sequence with %d elements", size);
	return NULL;
      }
      return (double *)PyArray_DATA(arrInput);
    }
  }
}


%typemap(in,fragment="getOrMakeArrayOfDoubles") double invector[ANY] (double temp[$1_dim0]) {
  $1 = getOrMakeArrayOfDoubles($input, $1_dim0, temp);
  if ($1 == NULL)
    return NULL;
  
}

%typemap(in,numinputs=0,noblock=1) double argoutvector[ANY] (PyObject *outArrObj), double &OUTPUT (PyObject *outArrObj) {
  int dims[1];
  dims[0] = 3;
  
  outArrObj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  
  if (outArrObj == NULL)
    return NULL;
  $1 = (double *)(PyArray_DATA((PyArrayObject *)outArrObj));
}

%typemap(argout,noblock=1) double argoutvector[ANY], double &argoutvector {
  $result = SWIG_Python_AppendOutput($result, outArrObj$argnum);
}

%init %{
  import_array();
%}
