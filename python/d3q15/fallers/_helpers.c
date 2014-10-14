#include "_helpers.h"

#include "../CDelta.h"
#include "../XArray.h"
#include "../utils.h"

#define NO_IMPORT_ARRAY
#include "../DQarrayobject.h"
#include <math.h>

/**
 * Update the positions of an array of fallers (periodic boundary
 * conditions).
 */
void FallerArray_move(PyObject *self, Lattice *lat) {
  PyObject *fallerData;
  double *fallers;
  int numFallers, numCols;
  int i;
  /* Start positions of particle data within each row of the table */
  int aStart, rStart, sStart, vStart, fStart;
  double hydroRadius, eta, noiseStDev;
  int size[DQ_d], tracks;
  
  /* get the Lattice size */
  size[0] = lat->nx;
  size[1] = lat->ny;
  size[2] = lat->nz;

  eta = lat->tau_s / 3.;
  
  /* abusing the name fallerData variable to store a ref to 
   * self.hydroRadius & .noiseStDev */
  fallerData = PyObject_GetAttrString(self, "hydroRadius");
  hydroRadius = PyFloat_AS_DOUBLE(fallerData);
  Py_DECREF(fallerData); /* discard the new ref gained above */

  fallerData = PyObject_GetAttrString(self, "noiseStDev");
  noiseStDev = PyFloat_AS_DOUBLE(fallerData);
  Py_DECREF(fallerData); /* discard the new ref gained above */

  fallerData = PyObject_GetAttrString(self, "tracks");
  tracks = (fallerData == Py_True);
  Py_DECREF(fallerData); /* discard the new ref gained above */

  /* now put fallerData to its proper use */
  fallerData = PyObject_GetAttrString(self, "data");
  /* noting we have a new ref *
   * and get the number of fallers */
  numFallers = PyArray_DIM((PyArrayObject *)fallerData, 0);
  /* and number of cols in array */
  numCols = PyArray_DIM((PyArrayObject *)fallerData, 1);

  fallers = (double *)PyArray_DATA( (PyArrayObject *)fallerData );
  rStart = XArray_getStart(self, 'r');
  if (tracks)
    sStart = XArray_getStart(self, 's');
  vStart = XArray_getStart(self, 'v');
  fStart = XArray_getStart(self, 'F');
  aStart = XArray_getStart(self, 'a');
  
#pragma omp parallel for schedule(guided)
  for (i=0; i<numFallers; ++i) {
    /* radius, position, velocity  & force on a particular particle */
    double a, *r, *s, *v, *F;
    double vinterp[DQ_d];
    int j;
    /* set up pointers to the correct parts of the array */
    r = fallers + numCols*i + rStart;
    if (tracks)
      s = fallers + numCols*i + sStart;
    v = fallers + numCols*i + vStart;
    F = fallers + numCols*i + fStart;
    a = fallers[numCols*i + aStart];
    
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
    

  } /* i */
  
  /* Drop the extra reference we picked up */
  Py_DECREF(fallerData);
}

void WalledFallerArray_move(PyObject *self, Lattice *lat) {
  PyObject *fallerData;
  double *fallers;
  int numFallers, numCols;
  int i;
  /* and their start positions within each row of the table */
  int aStart, rStart, sStart, vStart, fStart;
  double hydroRadius, eta, noiseStDev;
  int size[DQ_d], tracks;
  
  /* get the Lattice size */
  size[0] = lat->nx;
  size[1] = lat->ny;
  size[2] = lat->nz;

  eta = lat->tau_s / 3.;
  
  /* abusing the name fallerData variable to store a ref to 
   * self.hydroRadius & .noiseStDev */
  fallerData = PyObject_GetAttrString(self, "hydroRadius");
  hydroRadius = PyFloat_AS_DOUBLE(fallerData);
  Py_DECREF(fallerData); /* discard the new ref gained above */

  fallerData = PyObject_GetAttrString(self, "noiseStDev");
  noiseStDev = PyFloat_AS_DOUBLE(fallerData);
  Py_DECREF(fallerData); /* discard the new ref gained above */

  fallerData = PyObject_GetAttrString(self, "tracks");
  tracks = (fallerData == Py_True);
  Py_DECREF(fallerData); /* discard the new ref gained above */
  
  /* now put fallerData to its proper use */
  fallerData = PyObject_GetAttrString(self, "data");
  /* noting we have a new ref *
   * and get the number of fallers */
  numFallers = PyArray_DIM((PyArrayObject *)fallerData, 0);
  /* and number of cols in array */
  numCols = PyArray_DIM((PyArrayObject *)fallerData, 1);
  
  fallers = (double *)PyArray_DATA((PyArrayObject *)fallerData);
  rStart = XArray_getStart(self, 'r');
  if (tracks)
    sStart = XArray_getStart(self, 's');
  vStart = XArray_getStart(self, 'v');
  fStart = XArray_getStart(self, 'F');
  aStart = XArray_getStart(self, 'a');

#pragma omp parallel for schedule(guided)
  for (i=0; i<numFallers; ++i) {
    /* radius, position, velocity  & force on a particular particle */
    double a, *r, *s, *v, *F;
    double vinterp[DQ_d];
    int j;
    /* set up pointers to the correct parts of the array */
    r = fallers + numCols*i + rStart;
    if (tracks)
      s = fallers + numCols*i + sStart;
    v = fallers + numCols*i + vStart;
    F = fallers + numCols*i + fStart;
    a = fallers[numCols*i + aStart];
    
    /* if it's gone to the bottom already, leave it be */
    if (r[2] < 2.)
      continue;
	
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
    
    if (r[2] < 2.) {
      /* it has passed thru the bottom allowed position, 
       * so set its force to zero, so it won't affect things anymore */
      F[0] = 0.;
      F[1] = 0.;
      F[2] = 0.;
    }
    
  } /* i */
  
  /* Drop the extra reference we picked up */
  Py_DECREF(fallerData);
}

void PDFallerArray_addPotentialDipoles(PyObject *self, Lattice *lat) {
  PyObject *fallerData;
  double *fallers;
  int nx,ny,nz;
  int numFallers, numCols;

  int i,j,k,l,n,s; /* loop indices */
  int ind[3], size[3];
  double x[3], delta3d;
  /* radius, position of and force on a particular particle */
  double a, *r, *F;
  /* and their start positions within each row of the table */
  int aStart, rStart, fStart;
  double normF;
  double source_r[3], source_strength;
  double eta = lat->tau_s / 3.;
  
  /* get the PyObj of the data array */
  fallerData = PyObject_GetAttrString(self, "data");

  /* get pointers to the data sections */
  fallers = (double *)PyArray_DATA((PyArrayObject *)fallerData);

  /* get the Lattice size */
  size[0] = nx = lat->nx;
  size[1] = ny = lat->ny;
  size[2] = nz = lat->nz;
  
  /* and number of fallers */
  numFallers = PyArray_DIM((PyArrayObject *)fallerData, 0);
  /* and number of cols in array */
  numCols = PyArray_DIM((PyArrayObject *)fallerData, 1);
  
  rStart = XArray_getStart(self, 'r');
  fStart = XArray_getStart(self, 'F');
  aStart = XArray_getStart(self, 'a');

#define mod(a, n)  ((a)<1 ? (a)+n : ((a)>n ? (a)-n: (a)))
 
#pragma omp parallel for schedule(guided)
  for (n=0; n<numFallers; ++n) {
    /* set up pointers to the correct parts of the array */
    r = fallers + numCols*n + rStart;
    F = fallers + numCols*n + fStart;
    a = fallers[numCols*n + aStart];
    normF = sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);

    /* skip any unforced particle */
    if (normF == 0.)
      continue;

    for (s=-1; s<2; s+=2) {
      
      for  (l=0; l<3; l++) {
	source_r[l]  = r[l] + s*0.5 * a * F[l] / normF;
	if (source_r[l] < 0.5) {
	  source_r[l] += size[l];
	} else if (source_r[l] > size[l]+0.5) {
	  source_r[l] -= size[l];
	}
      }
      
      source_strength = s* -a * normF / (6. * eta);
      
      
      for (i= -2; i<2; ++i) {
	x[0] = ceil(source_r[0]+i);
	/* work out index of site given PBCs */
	ind[0] = mod(x[0], nx);
	
	for (j= -2; j<2; ++j) {
	  x[1] = ceil(source_r[1]+j);
	  ind[1] = mod(x[1], ny);
	  
	  for (k= -2; k<2; ++k) {
	    x[2] = ceil(source_r[2]+k);
	    ind[2] = mod(x[2], nz);
	    
	    /* evaluate the delta function */
	    delta3d = CDelta_delta(x[0] - source_r[0]) *
	      CDelta_delta(x[1] - source_r[1]) *
	      CDelta_delta(x[2] - source_r[2]);

	    for (l=0; l<DQ_q; ++l) 
	      DQ_f_get(lat, ind[0], ind[1], ind[2], l) += delta3d*source_strength/15.;
	    
	  } /* k */
	} /* j */
      } /* i */
      
    } /* sources, s */
    
  } /* fallers */
  
}
