#include "utils.h"
#define NO_IMPORT_ARRAY
#include "DQarrayobject.h"
#include "CDelta.h"

static int arrayOK(PyObject *o) {
  PyArrayObject *a;
  
  /* check it's a numpy array */
  if (!PyArray_Check(o)) { 
    PyErr_Format(PyExc_ValueError, "object must be a numpy array"); 
    return 0;
  }
  
  a = (PyArrayObject *)o;
  
  /* check array of doubles */
  if (PyArray_DTYPE(a)->type_num != NPY_DOUBLE) {
    PyErr_Format(PyExc_ValueError, "ArrayObject must be of type Float (a C double)");
    return 0;
  }
  
  /* check no. dimensions */
  if (PyArray_NDIM(a) != 2) {
    PyErr_Format(PyExc_ValueError, "ArrayObject must be 2-dimensional (is %d)", PyArray_NDIM(a));
    return 0;
  }
  
  /* check the 2nd dimension is 3D */
  if (PyArray_DIM(a, 1) != DQ_d) {
    PyErr_Format(PyExc_ValueError, "ArrayObject must be an array of %dD vectors (2nd array size is %" NPY_INTP_FMT ")", DQ_d, PyArray_DIM(a, 1));
    return 0;
  }
  
  return 1;
}

/* static double mod(double a, int n) { */
/*   if (a < 1.) { */
/*     return a+n; */
/*   } else if (a > n) { */
/*     return a-n; */
/*   } else { */
/*     return a; */
/*   } */
/* } */
void utils_addForceAtPosition_single(Lattice *lat,
				     double *F,
				     double *r) {
  int n[DQ_d];
  int indices[DQ_d][4];
  double deltas[DQ_d][4];
  double x, x0;
  int d;
  int i,j,k;

  n[DQ_X] = lat->nx;
  n[DQ_Y] = lat->ny;
  n[DQ_Z] = lat->nz;
  
  for (d=0; d< DQ_d; d++) {
    x0 = ceil(r[d]-2.);
    for (i=0; i<4; i++) {
      x = x0+i;
      indices[d][i] = (int)x;
      if (indices[d][i] < 1) {
	indices[d][i] += n[d];
      } else if (indices[d][i] > n[d]) {
	indices[d][i] -= n[d];
      }
      deltas[d][i] = CDelta_delta(r[d]-x);
    }
  }

  for (i=0; i<4; ++i) {
    for (j=0; j<4; ++j) {
      for (k=0; k<4; ++k) {
	
	/* add force contributions */
	for (d=0; d<DQ_d; ++d) {
	  DQ_force_get(lat,
		       indices[DQ_X][i],
		       indices[DQ_Y][j],
		       indices[DQ_Z][k],
		       d) += (deltas[DQ_X][i] *
			      deltas[DQ_Y][j] *
			      deltas[DQ_Z][k] *
			      F[d]);
	}
	
      }/* k */
    }/* j */
  }/* i */
  
}

void utils_addForceAtPosition_array(Lattice *lat,
				    PyObject *FObj,
				    PyObject *rObj) {
  PyArrayObject *FArrayObj, *rArrayObj;
  int num;
  /* Single position of and force */
  double *r, *F;
  int rNext, FNext;
  int n;
  
  if (!arrayOK(FObj) || !arrayOK(rObj))
    return;
  
  rArrayObj = (PyArrayObject *)rObj;
  FArrayObj = (PyArrayObject *)FObj;
  
  /* get number of pos/force pairs */
  num = PyArray_DIM(rArrayObj, 0);
  /* check forceArrayObj has the same number */
  if (PyArray_DIM(FArrayObj, 0) != num) {
    PyErr_Format(PyExc_ValueError, "FArrayObj and rArrayObj must have the same length");
    return;
  }
  
  r = (double *)PyArray_DATA(rArrayObj);
  F = (double *)PyArray_DATA(FArrayObj);
  rNext = PyArray_STRIDE(rArrayObj, 0)/sizeof(double);
  FNext = PyArray_STRIDE(FArrayObj, 0)/sizeof(double);
  
  for (n=0; n<num; ++n) {
    utils_addForceAtPosition_single(lat, F, r);
    
    r += rNext;
    F += FNext;
  }
}

void utils_interp_single(Lattice *lat, double *r, double *v) {
  int n[DQ_d];
  int indices[DQ_d][4];
  double deltas[DQ_d][4];
  double delta3d;
  double x, x0;
  int d;
  int i,j,k;

  n[DQ_X] = lat->nx;
  n[DQ_Y] = lat->ny;
  n[DQ_Z] = lat->nz;

  /* zero the data section */
  v[0] = v[1] = v[2] = 0.0;

  for (d=0; d< DQ_d; d++) {
    x0 = ceil(r[d]-2.);
    for (i=0; i<4; i++) {
      x = x0+i;
      indices[d][i] = (int)x;
      if (indices[d][i] < 1) {
	indices[d][i] += n[d];
      } else if (indices[d][i] > n[d]) {
	indices[d][i] -= n[d];
      }
      deltas[d][i] = CDelta_delta(r[d]-x);
    }
  }

  for (i=0; i<4; ++i) {
    for (j=0; j<4; ++j) {
      for (k=0; k<4; ++k) {
	/* evaluate the delta function */
	delta3d = deltas[DQ_X][i] * deltas[DQ_Y][j] * deltas[DQ_Z][k];
	
	for (d=0; d<3; ++d) {
	  v[d] += delta3d * DQ_u_get(lat,
				     indices[DQ_X][i],
				     indices[DQ_Y][j],
				     indices[DQ_Z][k],
				     d);
	}
      }
    }
  }
  
}
void utils_interp_array(Lattice *lat, PyObject *rObj, PyObject * vObj) {
  PyArrayObject *rArrayObj, *vArrayObj;
  int num, n;
  /* Single position of and force */
  double *r, *v;
  int rNext, vNext;
  
  if (!arrayOK(rObj) || !arrayOK(vObj))
    return;
  
  rArrayObj = (PyArrayObject *)rObj;
  vArrayObj = (PyArrayObject *)vObj;
  
  /* get number of pos/force pairs */
  num = PyArray_DIM(rArrayObj, 0);
  /* check forceArrayObj has the same number */
  if (PyArray_DIM(vArrayObj, 0) != num) {
    PyErr_Format(PyExc_ValueError, "rArrayObj and vArrayObj must have the same length");
    return;
  }

  r = (double *)PyArray_DATA(rArrayObj);
  v = (double *)PyArray_DATA(vArrayObj);
  rNext = PyArray_STRIDE(rArrayObj, 0)/sizeof(double);
  vNext = PyArray_STRIDE(vArrayObj, 0)/sizeof(double);
  
  for (n=0; n<num; ++n) {
    utils_interp_single(lat, r, v);
    r += rNext;
    v += vNext;
  } /* n */
  
}
