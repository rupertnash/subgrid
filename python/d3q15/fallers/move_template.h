/**
 * This is a pseudo template function with three named parameters that
 * must be #defined before including this file.
 *
 * FUNCTION_NAME - obvious
 *
 * BEGIN_ITER_CODE - code to be executed at the start of each
 * particle's move step.
 *
 * END_ITER_CODE - code to be executed at the end of each particle's
 * move step.
 */
void FUNCTION_NAME(PyObject *self, Lattice *lat) {
  PyObject *tmpObject;
  int numFallers;
  int i;
  /* XArray access context and member sub-arrays */
  XArrayCtx* ctx;
  XArrayMember aAr, rAr, sAr, vAr, fAr;

  double hydroRadius, eta, noiseStDev;
  int size[DQ_d], tracks;
  
  /* get the Lattice size */
  size[0] = lat->nx;
  size[1] = lat->ny;
  size[2] = lat->nz;

  eta = lat->tau_s / 3.;
  
  /* Use tmpObject to grab various attributes before conversion */
  tmpObject = PyObject_GetAttrString(self, "hydroRadius");
  hydroRadius = PyFloat_AS_DOUBLE(tmpObject);
  Py_DECREF(tmpObject); /* discard the new ref gained above */

  tmpObject = PyObject_GetAttrString(self, "noiseStDev");
  noiseStDev = PyFloat_AS_DOUBLE(tmpObject);
  Py_DECREF(tmpObject); /* discard the new ref gained above */

  tmpObject = PyObject_GetAttrString(self, "tracks");
  tracks = (tmpObject == Py_True);
  Py_DECREF(tmpObject); /* discard the new ref gained above */

  ctx = XArray_initCtx(self);
  XArray_getMember(ctx, 'r', &rAr);
  if (tracks)
    XArray_getMember(ctx, 's', &sAr);
  XArray_getMember(ctx, 'v', &vAr);
  XArray_getMember(ctx, 'F', &fAr);
  XArray_getMember(ctx, 'a', &aAr);
  
#pragma omp parallel for schedule(guided)
  for (i=0; i<ctx->nRows; ++i) {
    /* radius, position, velocity  & force on a particular particle */
    double a, *r, *s, *v, *F;
    double vinterp[DQ_d];
    int j;
    /* set up pointers to the correct parts of the array */
    r = XArray_getItem(&rAr, i);
    if (tracks)
      s = XArray_getItem(&sAr, i);
    v = XArray_getItem(&vAr, i);
    F = XArray_getItem(&fAr, i);
    a = *XArray_getItem(&aAr, i);
    
    BEGIN_ITER_CODE;
    
    /* interpolate the velocity to particle location */
    utils_interp_single(lat, r, vinterp);
    
    /* calculate velocity vector */
    for (j=0; j<DQ_d; ++j) {
      
      v[j] = vinterp[j] + /* advection */
	F[j] * (1./a - 1./hydroRadius) / (6. * M_PI * eta) + /* sedimentation*/
	noiseStDev*gasdev_get(lat->noise->gds); /* diffusion */

      /* Calculate the new position, assuming unit timestep */
      r[j] += v[j];
      if (tracks)
	s[j] += v[j];
      
      /*  Now check that it's not moved off the edge of the lattice */
      /*  and wrap it round if it has. */
      if (r[j] < 0.5) {
	r[j] += size[j];
      } else if (r[j] > size[j] + 0.5) {
	r[j] -= size[j];
      }
      
    }
    
    END_ITER_CODE;
  } /* i */
  
  /* Delete the array access ctx */
  XArray_delCtx(ctx);
}
