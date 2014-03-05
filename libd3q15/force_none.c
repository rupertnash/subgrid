#include "d3q15.h"

void force_none_init(Lattice *lat) {
  /* No setup needed, just zero the force field */
  for (int i=1; i<=lat->nx; i++)
    for (int j=1; j<=lat->ny; j++)
      for (int k=1; k<=lat->nz; k++)
	for (int d=0; d<DQ_d; d++)
	  DQ_force_get(lat, i,j,k, d) = 0.0;
}

void force_none_calc(Lattice *lat) {
  /* The distribution doesn't change, so do nothing. */
}
