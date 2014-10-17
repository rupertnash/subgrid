/* -*- mode: c; -*- */

%module d3q15
%{
#include "d3q15.h"
  /*numpy arrays */
#include "DQarrayobject.h"

#include "force_Py.c"
%}

%ignore Site;

%ignore mt_init;
%ignore mt_del;
%ignore mt_get;
%ignore mt_state::iv;

%ignore gasdev_init;
%ignore gasdev_del;
%ignore gasdev_get;
%ignore pgasdev_init;
%ignore pgasdev_del;
%ignore pgasdev_get;

%include ../../libd3q15/prng.h

%extend mt_state {
  int _bufsize() {
    return (2+DQ_PRNG_NTAB)*sizeof(long);
  }
  void _pickle(char* state) {
    //char* ans = (char*)malloc((2+DQ_PRNG_NTAB)*sizeof(long));
    long* buf = (long*)state;
    int i;
    for (i=0; i<DQ_PRNG_NTAB; i++)
      buf[i] = $self->iv[i];
    
    buf[i] = self->idum;
    buf[i+1] = self->iy;
  }
  void _unpickle(char* state) {
    long* buf = (long*)state;
    int i;
    for (i=0; i<DQ_PRNG_NTAB; i++)
      $self->iv[i] = buf[i];
    $self->idum = buf[i];
    $self->iy = buf[i+1];
  }
}
%extend gasdev_state {
  int _bufsize() {
    return sizeof(int) + sizeof(double) + mt_state__bufsize(NULL);
  }
  void _pickle(char* state) {
    *(int*)state = $self->iset;
    state += sizeof(int);
    *(double*)state = $self->gset;
    state += sizeof(double);
    mt_state__pickle(&$self->mt_st, state);
  }
  void _unpickle(char* state) {
    $self->iset = *(int*)state;
    state += sizeof(int);
    $self->gset = *(double*)state;
    state += sizeof(double);
    mt_state__unpickle(&$self->mt_st, state);
  }
}

%extend pgasdev_state {
   int _bufsize() {
     return sizeof(int) + $self->n * gasdev_state__bufsize(NULL);
   }
   void _pickle(char* state) {
     int i;
     *(int*)state = $self->n;
     state += sizeof(int);
     for (i = 0; i < $self->n; i++) {
       gasdev_state__pickle($self->ptrs+i, state);
       state += gasdev_state__bufsize($self->ptrs+i);
     }
   }
  void _unpickle(char* state) {
    int i;
    $self->n = *(int*)state;
    state += sizeof(int);
    for (i = 0; i < $self->n; i++) {
      gasdev_state__unpickle($self->ptrs+i, state);
      state += gasdev_state__bufsize($self->ptrs+i);
    }
  }
}
%ignore noise_init;
%ignore noise_set_temperature;
%ignore noise_add_to_modes;
%ignore noise_del;

%ignore NoiseConfig::rootT;
%ignore NoiseConfig::var;
%ignore NoiseConfig::_var;
%ignore NoiseConfig::temperature;
%ignore NoiseConfig::gds;

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
  int _bufsize() {
    return 12*sizeof(double) + pgasdev_state__bufsize(&$self->gds);
  }
  void _pickle(char* state) {
    int i;
    double* data = (double*)state;
    for (i=0; i<5; i++)
      *(data++) = $self->var[i];
    for (i=0; i<5; i++)
      *(data++) = $self->_var[i];
    *(data++) = $self->temperature;
    *(data++) = $self->rootT;
    
    pgasdev_state__pickle(&$self->gds, state + 12*sizeof(double));
  }
  void _unpickle(char* state) {
    int i;
    double* data = (double*)state;
    for (i=0; i<5; i++)
      $self->var[i] = *(data++);
    for (i=0; i<5; i++)
      $self->_var[i] = *(data++);
    $self->temperature = *(data++);
    $self->rootT = *(data++);
    
    pgasdev_state__unpickle(&$self->gds, state + 12*sizeof(double));
  }

  %pythoncode %{
    __swig_setmethods__["temperature"] = _d3q15.NoiseConfig_temperature_set
    __swig_getmethods__["temperature"] = _d3q15.NoiseConfig_temperature_get
    if _newclass:temperature = _swig_property(_d3q15.NoiseConfig_temperature_get, _d3q15.NoiseConfig_temperature_set)
    
    def _serialise(self):
        size = self._bufsize()
        state = size * ' '
        self._pickle(state)
        return state
    def _deserialise(self, state):
        size = self._bufsize()
        assert size == len(state), "Size of pickled state does not match required. Are you running on a different number of OpenMP threads?"
        self._unpickle(state)
  %}
}


%include Lattice.i
