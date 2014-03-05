/* -*- mode: c; -*- */

%include simpleNumpyTypemap.i

%module utils %{
#include "DQarrayobject.h"
#include "utils.h"
%}

%feature("autodoc", "1");
%pythoncode %{
import numpy
%}

%apply double invector[ANY] {double r[3]};
%apply double argoutvector[ANY] {double v[3]};

%rename(addForceAtPosition) utils_addForceAtPosition_single;
%rename(addForcesAtPositions) utils_addForceAtPosition_array;
%rename(interpOne) utils_interp_single;
%rename(interpArray) utils_interp_array;

%init %{
  import_array();
%}

%include utils.h

%pythoncode %{
    def interpArray(self, rArray, vArray=None):
        if vArray is None:
            vArray = numpy.zeros(rArray.shape)
        _utils.interpArray(self, rArray, vArray)
        return vArray
%}
