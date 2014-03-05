
void force_Py_calc(Lattice *lat) {
  /* The prototype of this function is fixed by the C library */

  PyObject *PyLat, *force_func_obj, *arglist;

  /* retrieve the PyObject of our Lattice */
  PyLat = (PyObject *)lat->config;
  
  /* from that, get the force_func_obj - note we have a new ref to it! */
  force_func_obj = PyObject_GetAttrString(PyLat, "__force_func_obj__");
  
  if (force_func_obj == NULL) {
    PyErr_SetString(PyExc_ValueError,
		 "Cannot retrieve Python function to calculate the force");
    return;
  }
  
  /* create the argument list, noting that we have a new ref to the tuple */
  arglist = Py_BuildValue("(O)", PyLat);
  
  /* call the function */
  PyEval_CallObject(force_func_obj, arglist);
  
  /* release our references */
  Py_DECREF(force_func_obj);
  Py_DECREF(arglist);
}
