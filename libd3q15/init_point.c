#include "d3q15.h"

void init_point (Lattice *lat, double rho, double u[DQ_d],
		 double rho_prime, int x[DQ_d]) {
  /* Initialise the distibution function arrays,
   * using u0 and v0 as the fluid velocities
   */
  double u_prime[DQ_d];
  
  init_uniform(lat, rho, u);
  
  /* increase rho at point specified by (x,y) */
  DQ_rho_get(lat, x[0],x[1],x[2]) = rho_prime;
  for (int d=0; d<DQ_d; d++)
    u_prime[d] = u[d] - 0.5*DQ_force_get(lat, x[0],x[1],x[2],d);
  
  calc_equil(rho_prime, u_prime, &DQ_f_get(lat, x[0],x[1],x[2],0));
  
  // done initialising
}
