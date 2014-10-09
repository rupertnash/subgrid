#include "d3q15.h"

void bc_pbc_init(Lattice *lat) {
  /* Do nothing (except check we're the right boundary condition)
   * but have to provide the function.
   */
  
}

void bc_pbc_update(Lattice *lat) {
  /* Updates the halo, giving periodic boundary conditions,
   * in all directions.
   */
  int i[DQ_d], n[DQ_d];
  
  n[0] = lat->nx; n[1] = lat->ny; n[2] = lat->nz;
  
  /* XY */
  for (i[0] = 0; i[0] <= n[0]+1; i[0]++) {
    for (i[1] = 0; i[1] <= n[0]+1; i[1]++) {
      i[2] = 0;
      bc_pbc_do_site(lat, i, n);
      i[2] = n[2] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
  /* XZ */
  for (i[0] = 0; i[0] <= n[0]+1; i[0]++) {
    for (i[2] = 1; i[2] <=  n[2]; i[2]++) {
      i[1] = 0;
      bc_pbc_do_site(lat, i, n);
      i[1] = n[1] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
  /* YZ */
  for (i[1] = 1; i[1] <= n[1]; i[1]++) {
    for (i[2] = 1; i[2] <= n[2]; i[2]++) {
      i[0] = 0;
      bc_pbc_do_site(lat, i, n);
      i[0] = n[0] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
    
  
}

#define isbulk(i) (i[0]>0 && i[0]<n[0]+1 && \
		   i[1]>0 && i[1]<n[1]+1 && \
		   i[2]>0 && i[2]<n[2]+1)

/* i = halo site */
void bc_pbc_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d]) {
  int d, p;
  int tgt[DQ_d], src[DQ_d];

  for (p=0; p< DQ_q; p++) {
    /* Work out which indices correspond to real f's that we got from
       the bulk. */
    int pPrime = lat->complement[p];
    for (d=0; d<DQ_d; d++) {
      src[d] = i[d]+ lat->ci[pPrime][d];
    }
    
    if isbulk(src) {
      /* We got f(i, p) from src - must send it to the other side of
       * the box. */
	
      for (d=0; d<DQ_d; d++) {
	tgt[d] = (src[d] - 1 + lat->ci[p][d]) % n[d] + 1;
      }
      
      DQ_f_get(lat, tgt[0],tgt[1],tgt[2], p) = DQ_f_get(lat, i[0],i[1],i[2], p);
    }
  }
}
#undef isbulk
