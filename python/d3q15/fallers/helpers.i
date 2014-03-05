/* -*- mode: c; -*- */

%module _helpers
%{
#include "../DQarrayobject.h"
#include "_helpers.h"
%}

%include _helpers.h

%init %{
  import_array();
%}
