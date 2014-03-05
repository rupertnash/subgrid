#include "d3q15.h"

void total_mass_and_momentum(Lattice *lat, double *mass, double mom[DQ_d]) {
  Site site;
  int i,j,k;
  int p;
  
  mass[0] = 0.;
  mom[DQ_X] = 0.;
  mom[DQ_Y] = 0.;
  mom[DQ_Z] = 0.;
  
  for (i=1; i<=lat->nx; i++) {
    for (j=1; j<=lat->ny; j++) {
      for (k=1; k<=lat->nz; k++) {
	set_site(lat, site, i,j,k);
	
	for (p=0; p<DQ_q; p++) {
	  mass[0] += site.f[p];
	  
	  mom[DQ_X] += lat->xi[p][DQ_X] * site.f[p];
	  mom[DQ_Y] += lat->xi[p][DQ_Y] * site.f[p];
	  mom[DQ_Z] += lat->xi[p][DQ_Z] * site.f[p];
	  
	}
	
      }
    }
  }
  
}
/*   double mom_plane_tot[DQ_d], mom_row_tot[DQ_d]; */
/*   double mass_plane_tot = 0.0, */
/*     mass_row_tot = 0.0; */
/*   double rho; */
/*   double rho_u[DQ_d]; */
/*   int i,j,k; */
/*   int d, p; */
/*   Site site; */
  
/*   *mass = 0.0; */
/*   for (d=0; d<DQ_d; d++) */
/*     mom[d] = 0.0; */
  
/*   for (i=1; i<=lat->nx; i++) { */
/*     mass_plane_tot = 0.0; */
/*     for (d=0; d<DQ_d; d++) */
/*       mom_plane_tot[d] = 0.0; */
    
/*     for (j=1; j<=lat->ny; j++) { */
/*       mass_row_tot = 0.0; */
/*       for (d=0; d<DQ_d; d++) */
/* 	mom_row_tot[d] = 0.0; */
      
/*       for (k=1; k<=lat->nz; k++) { */
	
/* 	set_site(lat, site, i,j,k); */

/* 	rho = 0.0; */
/* 	for (d=0; d<DQ_d; d++) */
/* 	  rho_u[d] = 0.; */
	
/* 	for (p=0; p<DQ_q; p++) { */
/* 	  rho += site.f[p] * lat->mm[DQ_rho][p]; */
/* 	  rho_u[DQ_X] += site.f[p] * lat->mm[DQ_momX][p]; */
/* 	  rho_u[DQ_Y] += site.f[p] * lat->mm[DQ_momY][p]; */
/* 	  rho_u[DQ_Z] += site.f[p] * lat->mm[DQ_momZ][p]; */
/* 	} */
	
/* 	mass_row_tot += rho; */
/* 	for (d=0; d<DQ_d; d++) */
/* 	  mom_row_tot[d] += rho_u[d] + 0.5 * rho * site.force[d]; */
/*       } */
      
/*       mass_plane_tot += mass_row_tot; */
/*       for (d=0; d<DQ_d; d++) */
/* 	mom_plane_tot[d] += mom_row_tot[d]; */
/*     } */
    
/*     *mass += mass_plane_tot; */
/*     for (d=0; d<DQ_d; d++) */
/*       mom[d] += mom_plane_tot[d]; */
/*   } */
/* } */
