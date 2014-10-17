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
#define FUNCTION_NAME FallerArray_move
#define BEGIN_ITER_CODE
#define END_ITER_CODE
#include "move_template.h"
#undef FUNCTION_NAME
#undef BEGIN_ITER_CODE
#undef END_ITER_CODE

#define FUNCTION_NAME WalledFallerArray_move
/* if it's gone to the bottom already, leave it be */
#define BEGIN_ITER_CODE if (r[2] < 2.) continue;
/* it has passed thru the bottom allowed position, 
 * so set its force to zero, so it won't affect things anymore */
#define END_ITER_CODE  if (r[2] < 2.) { F[0] = 0.; F[1] = 0.; F[2] = 0.; }
#include "move_template.h"
#undef FUNCTION_NAME
#undef BEGIN_ITER_CODE
#undef END_ITER_CODE

void PDFallerArray_addPotentialDipoles(PyObject *self, Lattice *lat) {
  int nx,ny,nz;
  int n, iPart;
  int size[3];
  /* start positions within each row of the table */
  int aStart, rStart, fStart;
  double eta = lat->tau_s / 3.;
  /* XArray access context and member sub-arrays */
  XArrayCtx* ctx;
  XArrayMember aAr, rAr, fAr;
  
  /* get the Lattice size */
  size[0] = nx = lat->nx;
  size[1] = ny = lat->ny;
  size[2] = nz = lat->nz;
    
  ctx = XArray_initCtx(self);
  XArray_getMember(ctx, 'r', &rAr);
  XArray_getMember(ctx, 'F', &fAr);
  XArray_getMember(ctx, 'a', &aAr);

#define mod(a, n)  ((a)<1 ? (a)+n : ((a)>n ? (a)-n: (a)))
 
#pragma omp parallel for schedule(guided)
  for (iPart=0; iPart<ctx->nRows; ++iPart) {
    int i,j,k,l,n,s; /* loop indices */
    int ind[3], size[3];
    double x[3], delta3d;
    /* radius, position of and force on a particular particle */
    double a, *r, *F;
    double normF;
    double source_r[3], source_strength;
    /* set up pointers to the correct parts of the array */
    r = XArray_getItem(&rAr, i);
    F = XArray_getItem(&fAr, i);
    a = *XArray_getItem(&aAr, i);

    normF = sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);

    /* skip any unforced particle */
    if (normF == 0.)
      continue;

    /* For each end */
    for (s=-1; s<2; s+=2) {
      /* For each dimension */
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
  XArray_delCtx(ctx);
}
