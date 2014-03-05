/* -*- mode: c; -*- */

%include simpleNumpyTypemap.i

%module(docstring="C implementation of Peskin's regularisation of the delta function.") CDelta %{
#include "CDelta.h"
%}

%feature("autodoc", "1");

%rename(delta) CDelta_delta;
%rename(delta3D) CDelta_3d_delta;

%apply double invector[ANY] {double x[ANY]}
%include CDelta.h
