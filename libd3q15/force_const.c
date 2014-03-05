#include "d3q15.h"

void force_const_init(Lattice *lat, double F[DQ_d]) {
  /* Set up the force array. It won't change during execution. */
  for (int i=1; i<=lat->nx; i++)
    for (int j=1; j<=lat->ny; j++)
      for (int k=1; k<=lat->nz; k++)
	for (int d=0; d<DQ_d; d++)
	  DQ_force_get(lat, i,j,k, d) = F[d];
}

void force_const_calc(Lattice *lat) {
  /* We don't actually need to do anything here,
   * though it mightn't be a bad idea to check that the force
   * hasn't bee altered by accident...
   */
}
