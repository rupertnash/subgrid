/* -*- mode: c; -*- */

%module d3q15
%{
#include "d3q15.h"
  /*numpy arrays */
#include "DQarrayobject.h"

#include "force_Py.c"
%}

%ignore Site;

%ignore noise_init;
%ignore noise_calc;
%ignore noise_set_temperature;
%ignore noise_add_to_modes;
%ignore noise_del;
%ignore gasdev;
%ignore ran1;

%ignore NoiseConfig::rootT;
%ignore NoiseConfig::var;
%ignore NoiseConfig::_var;
%ignore NoiseConfig::temperature;

%include ../../libd3q15/noise.h

%init %{
  /* To allow use of NumPy arrays */
  import_array();
%}

%extend NoiseConfig {
  void temperature_set(double t) {
    noise_set_temperature($self, t);
  }
  double temperature_get() {
    return $self->temperature;
  }

%pythoncode %{
    __swig_setmethods__["temperature"] = _d3q15.NoiseConfig_temperature_set
    __swig_getmethods__["temperature"] = _d3q15.NoiseConfig_temperature_get
    if _newclass:temperature = _swig_property(_d3q15.NoiseConfig_temperature_get, _d3q15.NoiseConfig_temperature_set)
%}
}
%include Lattice.i
