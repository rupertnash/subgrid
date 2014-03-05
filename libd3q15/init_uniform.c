#include "d3q15.h"

void init_uniform (Lattice *lat, double rho, double u[DQ_d]) {
  /* Initialise the distibution function arrays,
   * using u0 and v0 as the fluid velocities
   */
  double u_prime[DQ_d];
  for (int i=1; i<=lat->nx; i++) {
    for (int j=1; j<=lat->ny; j++) {
      for (int k=1; k<=lat->nz; k++) {
	/* we need to know the force and account for it to make sure
	 * that calc_hydro() correctly gives the initial conditions before
	 * any timesteps happen.
	 */ 
	DQ_rho_get(lat, i,j,k) = rho;
	
	for (int d=0; d<DQ_d; d++) {
	  DQ_u_get(lat, i,j,k,d) = u[d];
	  u_prime[d] = u[d] - 0.5*DQ_force_get(lat, i,j,k,d);
	}
	
	calc_equil(rho, u_prime, &DQ_f_get(lat, i,j,k,0));
      }
    }
  }
  // done initialising
}
