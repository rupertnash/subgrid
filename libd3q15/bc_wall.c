#include "d3q15.h"
#include <math.h>

/* static int complement[DQ_q]; */
static double wall_up_amp[DQ_d], wall_up_off[DQ_d],
  wall_lo_amp[DQ_d], wall_lo_off[DQ_d],
  wall_freq, wall_phase;

void bc_wall_init(Lattice *lat,
		  double v_up_amp[DQ_d], double v_up_off[DQ_d],
		  double v_lo_amp[DQ_d], double v_lo_off[DQ_d],
		  double shear_freq, double shear_phase) {
  /* Create the lookup table for which velocity is the opposite 
   * of each of the discrete velocities.
   */
  int d;
  
  for (d=0; d<DQ_d; d++) {
    wall_up_amp[d] = v_up_amp[d];
    wall_up_off[d] = v_up_off[d];
    wall_lo_amp[d] = v_lo_amp[d];
    wall_lo_off[d] = v_lo_off[d];
  }
  wall_freq = shear_freq;
  wall_phase = shear_phase;
  
}


#define isbulk(i) (i[0]>0 && i[0]<n[0]+1 && \
		   i[1]>0 && i[1]<n[1]+1 && \
		   i[2]>0 && i[2]<n[2]+1)

void bc_wall_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d], double v[DQ_d]) {
  /* We are at the top or bottom (i.e. the walls), so instead of 
   * doing pbc, we do bounce back, so copy the outward propagating 
   * distribution into the complementary velocity in the halo.
   */

  int d, p;
  int tgt[DQ_d];
  double cdotu, tmp, rho_site;
  
  for (p=1; p<DQ_q; p++) {
    
    /* For each velocity, work out the index of the node
     * that will be reached:
     *     tgt = i + xi[p]
     *
     * As a small optimization, don't consider the stationary
     * distribution, since it's not going anywhere.
     */
    for (d=0; d<DQ_d; d++) {
      tgt[d] = i[d]+  lat->xi[p][d];
    }
    
    if isbulk(tgt) {
      /* This velocity will propagate to a node in the bulk,
       * so do bounce back.
       */
      rho_site = DQ_rho_get(lat, tgt[0],tgt[1],tgt[2]);
      
      cdotu =  lat->xi[p][0]*v[0] +  lat->xi[p][1]*v[1] + lat->xi[p][2]*v[2];
      tmp = 2.0 * 3.0 * lat->w[p] * rho_site * cdotu ;
      
      DQ_f_get(lat, i[0],i[1],i[2], p) =
	DQ_f_get(lat, tgt[0],tgt[1],tgt[2], lat->complement[p]) + tmp;
      
    }
  }
  
  for (d=0; d<DQ_d; d++) 
    DQ_u_get(lat, i[0], i[1], i[2], d) = v[d];


}

void bc_wall_update(Lattice *lat) {
  /* Updates the halo, giving periodic boundary conditions,
   * in X & Y directions, walls in Z.
   */
  int i[DQ_d], n[DQ_d];
  double phase;
  double vlower[DQ_d], vupper[DQ_d];

  /* work out the possibly oscillating wall velocity at this time */
  if (wall_freq == 0.) {
    phase = 1.;
  } else {
    phase = sin(wall_freq* lat->time_step + wall_phase);
  }
  
  vlower[0] =  phase * wall_lo_amp[0] + wall_lo_off[0];
  vupper[0] =  phase * wall_up_amp[0] + wall_up_off[0];
  
  vlower[1] =  phase * wall_lo_amp[1] + wall_lo_off[1];
  vupper[1] =  phase * wall_up_amp[1] + wall_up_off[1];
  
  vlower[2] = vupper[2] = 0.;

  n[0] = lat->nx; n[1] = lat->ny; n[2] = lat->nz;
  
  /* XY  walls */
  for (i[0]=0; i[0]<=n[0]+1; i[0]++) {
    for (i[1]=0; i[1]<=n[1]+1; i[1]++) {
      i[2] = 0;
      bc_wall_do_site(lat, i, n, vlower);
      i[2] = n[2] + 1;
      bc_wall_do_site(lat, i, n, vupper);
    }
  }
  /* XZ periodic */
  for (i[0]=0; i[0]<=n[0]+1; i[0]++) {
    for (i[2]=1; i[2]<=n[2]; i[2]++) {
      i[1] = 0;
      bc_pbc_do_site(lat, i, n);
      i[1] = n[1] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
  /* YZ periodic */
  for (i[1]=1; i[1]<=n[1]; i[1]++) {
    for (i[2]=1; i[2]<=n[2]; i[2]++) {
      i[0] = 0;
      bc_pbc_do_site(lat, i, n);
      i[0] = n[0] + 1;
      bc_pbc_do_site(lat, i, n);
    }
  }
}
