#include "d3q15.h"
/* static int complement[DQ_q]; */

void bc_freeslip_init(Lattice *lat) {
}


#define isbulk(i) (i[0]>0 && i[0]<n[0]+1 && \
		   i[1]>0 && i[1]<n[1]+1 && \
		   i[2]>0 && i[2]<n[2]+1)


void bc_freeslip_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d]) {
  /* We are at the top or bottom (i.e. the walls), so instead of 
   * doing pbc, we do bounce back, so copy the outward propagating 
   * distribution into the complementary velocity in the halo.
   */
  int d, p;
  int tgt[DQ_d], ori[DQ_d];

  for (p=1; p< DQ_q; p++) {
    /* For each velocity, work out the index of the node
     * that will be reached:
     *     tgt = i + xi[p]
     *
     * As a small optimization, don't consider the stationary
     * distribution, since it's not going anywhere.
     */
    for (d=0; d<DQ_d; d++) {
      tgt[d] = i[d]+ lat->xi[p][d];
    }
    
    if isbulk(tgt) {
      /* This velocity will propagate to a node in the bulk,
       * work out where the originating node is
       */
      /* Note this only works for z walls */
      ori[0] = i[0] < 1 ? n[0] : (i[0] > n[0] ? 1 : i[0]);
      ori[1] = i[1] < 1 ? n[1] : (i[1] > n[1] ? 1 : i[1]);
      ori[2] = tgt[2];
      
      /* Work out which of velocities to use */
      /* because of ordering, can use (p-1) XOR 1, +1 */
      DQ_f_get(lat, i[0],i[1],i[2], p) = 
	DQ_f_get(lat, ori[0],ori[1],ori[2], ((p-1) ^ 1) + 1);
    } 
  }
}

void bc_freeslip_update(Lattice *lat) {
  /* Updates the halo, giving periodic boundary conditions,
   * in X & Y directions, walls in Z.
   */
  int i[DQ_d], n[DQ_d];

  n[0] = lat->nx; n[1] = lat->ny; n[2] = lat->nz;
  
  /* XY  walls */
  for (i[DQ_X]=0; i[DQ_X]<=n[DQ_X]+1; i[DQ_X]++) {
    for (i[DQ_Y]=0; i[DQ_Y]<=n[DQ_Y]+1; i[DQ_Y]++) {
      i[DQ_Z] = 0;
      bc_noslip_do_site(lat, i, n);
      i[DQ_Z] = n[DQ_Z] + 1;
      bc_freeslip_do_site(lat, i, n);
    }
  }
  /* XZ periodic */
  for (i[DQ_X]=0; i[DQ_X]<=n[DQ_X]+1; i[DQ_X]++) {
    for (i[DQ_Z]=1; i[DQ_Z]<=n[DQ_Z]; i[DQ_Z]++) {
      i[DQ_Y] = 0;
      bc_pbc_do_site(lat, i, n);
      i[DQ_Y] = n[DQ_Y] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
  /* YZ periodic */
  for (i[DQ_Y]=1; i[DQ_Y]<=n[DQ_Y]; i[DQ_Y]++) {
    for (i[DQ_Z]=1; i[DQ_Z]<=n[2]; i[DQ_Z]++) {
      i[DQ_X] = 0;
      bc_pbc_do_site(lat, i, n);
      i[DQ_X] = n[DQ_X] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
}
