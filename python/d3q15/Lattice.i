/* -*- mode: c; -*- */
%ignore Lattice::f_current_ptr;
%ignore Lattice::f_new_ptr;
%ignore Lattice::rho_ptr;
%ignore Lattice::force_ptr;
%ignore Lattice::u_ptr;
%ignore Lattice::w;
%ignore Lattice::xi;
%ignore Lattice::ci;
%ignore Lattice::Q;
%ignore Lattice::complement;
%ignore Lattice::mm;
%ignore Lattice::mmi;
%ignore Lattice::norms;

%ignore Lattice::bc_func;
%ignore Lattice::force_func;

%include ../../libd3q15/Lattice.h

/* SWIG macro to add an exception check to wrapper code.
 * Returning NULL is the standard indication to the Python 
 * API that something went wrong. Success of a void function 
 * will return a (new) reference to Py_None to the API. */

%define EXC_CHECK(method)
%exception method {
  $action
    if (PyErr_Occurred()) {
      return NULL;
      
    }
}
%enddef

/* any function which returns void and uses (even indirectly)
 * any Python API functions must be listed here to have 
 * exception checking done.
 */
EXC_CHECK(initForceC)
EXC_CHECK(setForceC)
EXC_CHECK(setForcePy)
EXC_CHECK(__initBC_periodic)
EXC_CHECK(__initBC_noslip)
EXC_CHECK(__initBC_freeslip)
EXC_CHECK(__initBC_wall)
EXC_CHECK(step)

EXC_CHECK(f_set)
EXC_CHECK(rho_set)
EXC_CHECK(u_set)
EXC_CHECK(force_set)

%extend Lattice {
  // Constructor
  Lattice(int nx, int ny, int nz, double tau_s, double tau_b) {
    Lattice *lat = d3q15_init(nx, ny, nz, tau_s, tau_b);
    return lat;
  }
  
  void __updateBackPointer__(PyObject *PyLat) {
    /* This function sets the config member of the C struct to point
     * to the PyObject that wraps the C data structure.
     * Slightly ugly to have the same object as self and as PyLat,
     * but needed to have SWIG do the conversion for us.
     */
    $self->config = (void *)PyLat;
  }
  
  // Destructor
  ~Lattice() {
    d3q15_destroy($self);
  }
  
  void updateHydroVars() {
    calc_hydro($self);
  }

  /************************************************************
   * "methods" to register the forcing function to the C
   * library functions
   ************************************************************/
  void initForceC(char *name, ...) {
    /* Initialises the force field and sets the forcing to
     * a C function, whose name is one of the options below.
     */
    va_list arg_p;
    va_start(arg_p, name);
    
    if (!strcmp(name, "none") ){
      force_none_init($self);
      $self->force_func = force_none_calc;
      
    } else if (!strcmp(name, "const")) {
      double F[DQ_d];
      F[0] = va_arg(arg_p, double);
      F[1] = va_arg(arg_p, double);
      F[2] = va_arg(arg_p, double);
      force_const_init($self, F);
      $self->force_func = force_const_calc;

    } else {
      PyErr_Format(PyExc_ValueError, "unknown forcing type %s", name);
    }
    
    va_end(arg_p);
  }
  
  void setForceC(char *name) {
    /* Just sets the forcing to a C function,
     * whose name is one of the options below.
     * DOES NOT INITIALISE!
     */
    if (!strcmp(name, "none") ){
      $self->force_func = force_none_calc;
      
    } else if (!strcmp(name, "const")) {
      $self->force_func = force_const_calc;
    } else {
      PyErr_Format(PyExc_ValueError, "unknown forcing type %s", name);
    }
    
  }
  
  void setForcePy(PyObject *PyLat, PyObject *calc_func) {
    /* This tells the C code to use the Python fucntion given to
     * update the force field when it needs it. Is wrapped by
     * force_Py_calc(Lattice *) defined in force_Py.c
     *
     * Note that the calling Python script is expected to
     * do the initialisation.
     *
     * IMPORTANT: Note the repetition of the lattice object in the 
     * argument list (first as Lattice *self, second as PyObject *PyLat)
     * This is to get SWIG to do pointer conversion for us.
     * Also, this method is overriden and wrapped in __init__.py to
     * make the interface cleaner.
     */
    
    /* Check that it's callable, if not, raise ValueError */  
    if (!PyCallable_Check(calc_func)) {
      PyErr_Format(PyExc_ValueError, "function supplied is not callable!");
    }
    
    /* Store the force function object in the Python object */
    PyObject_SetAttrString(PyLat, "__force_func_obj__", calc_func);
    
    /* Set the C function that updates the force to be the one in this file */
    $self->force_func = force_Py_calc;

  }
  
  %pythoncode %{
    def setForcePy(self, forceFunc):
        """Make the Python callable forceFunc be called at the appropriate
        moment by the C level code with the (Python level) Lattice object 
        as its only argument.
        
        This wraps the lower level function to give a cleaner interface."""
      
        return _d3q15.Lattice_setForcePy(self, self, forceFunc)
  %}

  void __initBC_periodic() {
    bc_pbc_init($self);
    $self->bc_func = bc_pbc_update;
  }
  
  void __initBC_noslip() {
    bc_noslip_init($self);
    $self->bc_func = bc_noslip_update;
  }
  
  void __initBC_freeslip() {
    bc_freeslip_init($self);
    $self->bc_func = bc_freeslip_update;
  }
  
  void __initBC_wall(double v_up_amp_x, double v_up_amp_y,
		     double v_up_off_x, double v_up_off_y,
		     double v_lo_amp_x, double v_lo_amp_y,
		     double v_lo_off_x, double v_lo_off_y,
		     double freq, double phase) {
    double v_up_amp[DQ_d];
    double v_up_off[DQ_d];
    double v_lo_amp[DQ_d];
    double v_lo_off[DQ_d];

    v_up_amp[0] = v_up_amp_x; v_up_amp[1] = v_up_amp_y; v_up_amp[2] = 0;
    v_up_off[0] = v_up_off_x; v_up_off[1] = v_up_off_y; v_up_off[2] = 0;
    v_lo_amp[0] = v_lo_amp_x; v_lo_amp[1] = v_lo_amp_y; v_lo_amp[2] = 0;
    v_lo_off[0] = v_lo_off_x; v_lo_off[1] = v_lo_off_y; v_lo_off[2] = 0;
    
    bc_wall_init($self, 
		 v_up_amp, v_up_off,
		 v_lo_amp, v_lo_off,
		 freq, phase);
    $self->bc_func = bc_wall_update;
  }

  %pythoncode %{
    def initBoundaryC(self, name, *args):
        """Initialize the boundary conditions to a C function specified by
        the string 'name'. Current options are:
        
        periodic -- full 3D PBCs
                    no further arguments
        
        noslip -- stationary solid walls in the z-direction, PBC in x- and y-
                  no further arguments
        
        wall -- tangentially moving, possibly oscillating solid walls in z, 
                PBC in x & y
                Arguments:
                v_up_amp: vector amplitude of oscillation of upper wall velocity
                v_up_off: vector constant part of upper wall velocity
                v_lo_amp: vector amplitude of oscillation of lower wall velocity
                v_lo_off: vector constant part of lower wall velocity
                freq -- angular frequency of oscillation
                phase -- phase shift in radians
                
                The actual velocity of the walls are calculated as:
                vlower =   v_lo_amp*sin(freq*time + phase) + wall_lo_off
                and similarly for vupper.
 
        freeslip -- free-slip walls in the z-direction, PBC in x- and y-
            no further arguments
        **** THIS LAST HAS NOT BEEN TESTED ****

        """
        if name == 'periodic':
            self.__initBC_periodic()
        elif name == 'noslip':
            self.__initBC_noslip()
        elif name == 'freeslip':
            self.__initBC_freeslip()
        elif name == 'wall':
            v_up_amp, v_up_off, v_lo_amp, v_lo_off, freq, phase = args
            self.__initBC_wall(v_up_amp[0], v_up_amp[1],
                               v_up_off[0], v_up_off[1],
                               v_lo_amp[0], v_lo_amp[1],
                               v_lo_off[0], v_lo_off[1],
                               freq, phase)
        else:
            raise ValueError('Unknown boundary type')
        
        return
  %}
  
  /* Given the hydrodynamic fields, initialises the distribution
   * function array to their equilibrium values.
   * Might not correctly deal with forcing however... need to check that!
   */
  void initFromHydroVars() {
    int i, j, k;
    for (i=1; i<=$self->nx; i++)
      for (j=1; j<=$self->ny; j++)
	for (k=1; k<=$self->nz; k++)
	  calc_equil(DQ_rho_get($self,i,j,k), &DQ_u_get($self, i,j,k,0), &DQ_f_get($self, i,j,k,0));
  }

  /* pseduo method wrapper */
  void step() {
    d3q15_step($self);
  }
  
  /* methods to get the sizes of the arrays */
  PyObject *scalarFieldSize() {
    return Py_BuildValue("iii",
			 $self->_x_ar_size, $self->_y_ar_size, $self->_z_ar_size);
  }
  PyObject *vectorFieldSize() {
    return Py_BuildValue("iiii",
			 $self->_x_ar_size, $self->_y_ar_size, $self->_z_ar_size,
			 DQ_d);
  }
  PyObject *fFieldSize() {
    return Py_BuildValue("iiii",
			 $self->_x_ar_size, $self->_y_ar_size, $self->_z_ar_size,
			 DQ_q);
  }
  
  const double totalMass;
  PyObject *totalMomentum;

  PyObject *f_current_get() {
    /* Return a NumPy array Object whose data section points
     * to the Lattice's distribution data.
     */
    PyObject *resultobj = NULL;
    int n_dims = DQ_d + 1;
    npy_intp dims[n_dims];
    
    dims[0] = $self->nx + 2;
    dims[1] = $self->ny + 2;
    dims[2] = $self->nz + 2;
    dims[3] = DQ_q;
    
    resultobj = PyArray_SimpleNewFromData(n_dims,
					  dims,
					  NPY_DOUBLE,
					  (void *)$self->f_current_ptr);
    
    return resultobj;
  }
  
  void f_current_set(PyObject *obj) {
    PyArrayObject *array = NULL;
    double *data = NULL;
    int i;
    
    /* check it's a numpy array */  
    if (!PyArray_Check(obj)) {
      PyErr_Format(PyExc_ValueError, "must set to a numpy array");
      return;
    }
    
    array = (PyArrayObject *)obj;
    /* check it's an array of doubles */
    if (PyArray_DTYPE(array)->type_num != NPY_DOUBLE) {
      PyErr_Format(PyExc_ValueError, "array must be of type Float (a C double)");
      return;
    }
    
    /* check it has the right no. dimensions */
    if (PyArray_NDIM(array) != DQ_d+1) {
      PyErr_Format(PyExc_ValueError, "array must be %d-dimensional (is %d)", DQ_d+1, PyArray_NDIM(array));
      return;
    }
    
    /* check the dimensions match */
    npy_intp *dims = PyArray_DIMS(array);
    if (dims[0] != $self->_x_ar_size || 
        dims[1] != $self->_y_ar_size ||
        dims[2] != $self->_z_ar_size ||
        dims[3] != DQ_q) {
      PyErr_Format(PyExc_ValueError,
		   "array must sized %d x %d x %d x %d (is %" NPY_INTP_FMT 
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT ")",
		   $self->_x_ar_size,$self->_y_ar_size,$self->_z_ar_size,DQ_q,
		   dims[0], dims[1],
		   dims[2], dims[3]);
      return;
    }
    
    /* copy */
    data = (double *)PyArray_DATA(array);
    for (i=0; i<($self->nx+2)*($self->ny+2)*($self->nz+2)*DQ_q; i++)
      $self->f_current_ptr[i] = data[i];
    
  }
  
  PyObject *f_new_get() {
    /* Return a NumPy array Object whose data section points
     * to the Lattice's distribution data.
     */
    PyObject *resultobj = NULL;
    int n_dims = DQ_d + 1;
    npy_intp dims[n_dims];
    
    dims[0] = $self->nx + 2;
    dims[1] = $self->ny + 2;
    dims[2] = $self->nz + 2;
    dims[3] = DQ_q;
    
    resultobj = PyArray_SimpleNewFromData(n_dims,
					  dims,
					  NPY_DOUBLE,
					  (void *)$self->f_new_ptr);
    
    return resultobj;
  }
  
  void f_new_set(PyObject *obj) {
    PyArrayObject *array = NULL;
    double *data = NULL;
    int i;
    
    /* check it's a numpy array */  
    if (!PyArray_Check(obj)) {
      PyErr_Format(PyExc_ValueError, "must set to a numpy array");
      return;
    }
    
    array = (PyArrayObject *)obj;
    /* check it's an array of doubles */
    if (PyArray_DTYPE(array)->type_num != NPY_DOUBLE) {
      PyErr_Format(PyExc_ValueError, "array must be of type Float (a C double)");
      return;
    }
    
    /* check it has the right no. dimensions */
    if (PyArray_NDIM(array) != DQ_d+1) {
      PyErr_Format(PyExc_ValueError, "array must be %d-dimensional (is %d)", DQ_d+1, PyArray_NDIM(array));
      return;
    }
    
    /* check the dimensions match */
    npy_intp *dims = PyArray_DIMS(array);
    if (dims[0] != $self->_x_ar_size || 
        dims[1] != $self->_y_ar_size ||
        dims[2] != $self->_z_ar_size ||
        dims[3] != DQ_q) {
      PyErr_Format(PyExc_ValueError,
		   "array must sized %d x %d x %d x %d (is %" NPY_INTP_FMT 
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT ")",
		   $self->_x_ar_size,$self->_y_ar_size,$self->_z_ar_size,DQ_q,
		   dims[0], dims[1],
		   dims[2], dims[3]);
      return;
    }
    
    /* copy */
    data = (double *)PyArray_DATA(array);
    for (i=0; i<($self->nx+2)*($self->ny+2)*($self->nz+2)*DQ_q; i++)
      $self->f_new_ptr[i] = data[i];
    
  }
  
  /******* RHO ********/
  PyObject *rho_get() {
    /* Return a NumPy array Object whose data section points
     * to the Lattice's density
     */
    PyObject *resultobj = NULL;
    int n_dims = DQ_d;
    npy_intp dims[n_dims];
    
    dims[0] = $self->nx + 2;
    dims[1] = $self->ny + 2;
    dims[2] = $self->nz + 2;
    
    resultobj = PyArray_SimpleNewFromData(n_dims,
					  dims,
					  NPY_DOUBLE,
					  (void *)$self->rho_ptr);
    
    return resultobj;
  }
  
  void rho_set(PyObject *obj) {
    PyArrayObject *array = NULL;
    double *data = NULL;
    int i;
    
    /* check it's a numpy array */  
    if (!PyArray_Check(obj)) {
      PyErr_Format(PyExc_ValueError, "must set to a numpy array");
      return;
    }
    
    array = (PyArrayObject *)obj;
    /* check it's an array of doubles */
    if (PyArray_DTYPE(array)->type_num != NPY_DOUBLE) {
      PyErr_Format(PyExc_ValueError, "array must be of type Float (a C double)");
      return;
    }
    
    /* check it has the right no. dimensions */
    if (PyArray_NDIM(array) != DQ_d) {
      PyErr_Format(PyExc_ValueError, "array must be %d-dimensional (is %d)",
		   DQ_d, PyArray_NDIM(array));
      return;
    }
    
    /* check the dimensions match */
    npy_intp *dims = PyArray_DIMS(array);
    if (dims[0] != $self->_x_ar_size || 
        dims[1] != $self->_y_ar_size ||
        dims[2] != $self->_z_ar_size) {


      PyErr_Format(PyExc_ValueError, "array must sized %d x %d x %d"
		   "(is %" NPY_INTP_FMT 
		   " x %" NPY_INTP_FMT 
		   " x %" NPY_INTP_FMT ")",
		   $self->_x_ar_size, $self->_y_ar_size, $self->_z_ar_size,
  		   dims[0], dims[1],
		   dims[2]);
      return;
    }
    
    /* copy */
    data = (double *)PyArray_DATA(array);
    for (i=0; i<($self->nx+2)*($self->ny+2)*($self->nz+2); i++)
      $self->rho_ptr[i] = data[i];
    
  }
  
  /************** U **********/
  PyObject *u_get() {
    /* Return a NumPy array Object whose data section points
     * to the Lattice's velocity field
     */
    PyObject *resultobj = NULL;
    int n_dims = DQ_d + 1;
    npy_intp dims[n_dims];
    
    dims[0] = $self->nx + 2;
    dims[1] = $self->ny + 2;
    dims[2] = $self->nz + 2;
    dims[3] = DQ_d;

    resultobj = PyArray_SimpleNewFromData(n_dims,
					  dims,
					  NPY_DOUBLE,
					  (void *)$self->u_ptr);
    
    return resultobj;
  }
  
  void u_set(PyObject *obj) {
    PyArrayObject *array = NULL;
    double *data = NULL;
    int i;
    
    /* check it's a numpy array */  
    if (!PyArray_Check(obj)) {
      PyErr_Format(PyExc_ValueError, "must set to a numpy array");
      return;
    }
    
    array = (PyArrayObject *)obj;
    /* check it's an array of doubles */
    if (PyArray_DTYPE(array)->type_num != NPY_DOUBLE) {
      PyErr_Format(PyExc_ValueError, "array must be of type Float (a C double)");
      return;
    }
    
    /* check it has the right no. dimensions */
    if (PyArray_NDIM(array) != DQ_d+1) {
      PyErr_Format(PyExc_ValueError, "array must be %d-dimensional (is %d)",
		   DQ_d+1, PyArray_NDIM(array));
      return;
    }
    
    /* check the dimensions match */
    npy_intp *dims = PyArray_DIMS(array);
    if (dims[0] != $self->_x_ar_size || 
	dims[1] != $self->_y_ar_size ||
	dims[2] != $self->_z_ar_size ||
	dims[3] != DQ_d) {
      PyErr_Format(PyExc_ValueError,
		   "array must sized %d x %d x %d x %d (is %" NPY_INTP_FMT 
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT ")",
		   $self->_x_ar_size, $self->_y_ar_size, $self->_z_ar_size, DQ_d,
		   dims[0], dims[1],
		   dims[2], dims[3]);
      return;
    }
    
    /* copy */
    data = (double *)PyArray_DATA(array);
    for (i=0; i<($self->nx+2)*($self->ny+2)*($self->nz+2)*DQ_d; i++)
      $self->u_ptr[i] = data[i];
    
  }
  
  /************** FORCE **********/
  PyObject *force_get() {
    /* Return a NumPy array Object whose data section points
     * to the Lattice's force field.
     */
    PyObject *resultobj = NULL;
    int n_dims = DQ_d + 1;
    npy_intp dims[n_dims];
    
    dims[0] = $self->nx + 2;
    dims[1] = $self->ny + 2;
    dims[2] = $self->nz + 2;
    dims[3] = DQ_d;
    
    resultobj = PyArray_SimpleNewFromData(n_dims,
					  dims,
					  NPY_DOUBLE,
					  (void *)$self->force_ptr);
    
    return resultobj;
  }
  
  void force_set(PyObject *obj) {
    PyArrayObject *array = NULL;
    double *data = NULL;
    int i;
    
    /* check it's a numpy array */  
    if (!PyArray_Check(obj)) {
      PyErr_Format(PyExc_ValueError, "must set to a numpy array");
      return;
    }
    
    array = (PyArrayObject *)obj;
    /* check it's an array of doubles */
    if (PyArray_DTYPE(array)->type_num != NPY_DOUBLE) {
      PyErr_Format(PyExc_ValueError, "array must be of type Float (a C double)");
      return;
  }
    
    /* check it has the right no. dimensions */
    if (PyArray_NDIM(array) != DQ_d+1) {
      PyErr_Format(PyExc_ValueError, "array must be %d-dimensional (is %d)",
		   DQ_d+1, PyArray_NDIM(array));
      return;
    }
    
    /* check the dimensions match */
    npy_intp *dims = PyArray_DIMS(array);
    if (dims[0] != $self->_x_ar_size || 
	dims[1] != $self->_y_ar_size ||
	dims[2] != $self->_z_ar_size ||
	dims[3] != DQ_d) {
      PyErr_Format(PyExc_ValueError,
		   "array must sized %d x %d x %d x %d (is %" NPY_INTP_FMT 
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT
		   " x %" NPY_INTP_FMT ")",
		   $self->_x_ar_size,$self->_y_ar_size,$self->_z_ar_size,DQ_d,
		   dims[0], dims[1],
		   dims[2], dims[3]);
      return;
    }
    
    /* copy */
    data = (double *)PyArray_DATA(array);
    for (i=0; i<($self->nx+2)*($self->ny+2)*($self->nz+2)*DQ_d; i++)
      $self->force_ptr[i] = data[i];
    
  }
  
/* This python code, inserted into the shadow class, uses the getters & 
 * setters to implement properties for the distributions, rho, u & force
 */
%pythoncode %{
    __swig_setmethods__["f_current"] = _d3q15.Lattice_f_current_set
    __swig_getmethods__["f_current"] = _d3q15.Lattice_f_current_get
    if _newclass:f_current = _swig_property(_d3q15.Lattice_f_current_get, _d3q15.Lattice_f_current_set)
    
    __swig_setmethods__["rho"] = _d3q15.Lattice_rho_set
    __swig_getmethods__["rho"] = _d3q15.Lattice_rho_get
    if _newclass:rho = _swig_property(_d3q15.Lattice_rho_get, _d3q15.Lattice_rho_set)
    
    __swig_setmethods__["u"] = _d3q15.Lattice_u_set
    __swig_getmethods__["u"] = _d3q15.Lattice_u_get
    if _newclass:u = _swig_property(_d3q15.Lattice_u_get, _d3q15.Lattice_u_set)
    
    __swig_setmethods__["force"] = _d3q15.Lattice_force_set
    __swig_getmethods__["force"] = _d3q15.Lattice_force_get
    if _newclass:force = _swig_property(_d3q15.Lattice_force_get, _d3q15.Lattice_force_set)
%}
  
%pythoncode %{
    def __reduce__(self):
        """For pickling. This returns a tuple (callable, args, state)
        defined in section 13.1.5.2 of the Python Library Reference."""
        initArgs = (self.nx, self.ny, self.nz, self.tau_s, self.tau_b)
        return (self.__class__, initArgs, self.__getstate__())
    
    def __getstate__(self):
        """Enables pickling - creates the actual data to be
        serialized."""
        picdict = {}
        picdict['fields'] = {'f_current': self.f_current,
			     'f_new': self.f_new,
                             'force': self.force,
                             'u': self.u,
                             'rho': self.rho}
        picdict['time_step'] = self.time_step
        picdict['temperature'] = self.noise.temperature
        picdict['noiseSeed'] = self.noise.seed
        return picdict
    
    def __setstate__(self, picdict):
        """Restores from the pickled data."""
        self.time_step = picdict['time_step']
        self.noise.temperature = picdict['temperature']
        self.noise.seed = picdict['noiseSeed']
        self.f_current = picdict['fields']['f_current']
        self.f_new = picdict['fields']['f_new']
        self.force = picdict['fields']['force']
        self.u = picdict['fields']['u']
        self.rho = picdict['fields']['rho']
        return

%}

}

%{

double Lattice_totalMass_get(Lattice *self) {
  double mass, mom[DQ_d];
  total_mass_and_momentum(self, &mass, mom);
  return mass;
}

PyObject *Lattice_totalMomentum_get(Lattice *self) {
  int i;
  npy_intp dims = DQ_d;
  double mass, mom[DQ_d], *ans_data;
  PyObject *ans = PyArray_SimpleNew(1, &dims, NPY_DOUBLE);
  PyArrayObject* ansArray = (PyArrayObject *)ans;
  
  total_mass_and_momentum(self, &mass, mom);
  ans_data = (double *)PyArray_DATA(ansArray);
  for (i=0; i<DQ_d; i++)
    ans_data[i] = mom[i];
  
  return ans;
}
 
void Lattice_totalMomentum_set(Lattice *self, PyObject *arg) {
  
}


%}
			
