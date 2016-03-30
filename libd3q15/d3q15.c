/* 3 dimensional, 15 velocity Lattice Boltzmann code.
 */
#include <stdlib.h>
#include <math.h>
#include "d3q15.h"

/************************************************************/

static const int q = DQ_q;
static const int d = DQ_d;

Lattice *d3q15_init(int nx, int ny, int nz, double tau_s, double tau_b) {
  Lattice *lat = (Lattice *)malloc(sizeof(Lattice));
  
  /* Here we include the out put of eigenvectors.py */
#define DQ_init
#include "eigenvectors.h"
#undef DQ_init

  lat->nx = nx;
  lat->ny = ny;
  lat->nz = nz;
  lat->_x_ar_size = nx + 2;
  lat->_y_ar_size = ny + 2;
  lat->_z_ar_size = nz + 2;
  
  lat->strides[DQ_X] = (ny+2) * (nz+2);
  lat->strides[DQ_Y] = (nz+2);
  lat->strides[DQ_Z] = 1;
  
  /* Allocate space for our arrays */
  lat->f_ptr = (double *)malloc( (nx+2) * (ny+2) * (nz+2) * DQ_q *
				 sizeof(double) );
  lat->rho_ptr = (double *)malloc( (nx+2) * (ny+2) * (nz+2) *
				   sizeof(double) );
  lat->u_ptr = (double *)malloc( (nx+2) * (ny+2) * (nz+2) * DQ_d *
				 sizeof(double) );
  
  lat->force_ptr = (double *)malloc( (nx+2) * (ny+2) * (nz+2) * DQ_d *
				     sizeof(double) );
  lat->tau_b = tau_b;
  lat->tau_s = tau_s;
  
  lat->time_step = 0;
  
  /* Set our function pointers to NULL so we can know if they 
   * are initialised or not. */
  lat->force_func = NULL;
  lat->bc_func = NULL;

#ifdef DQ_NOISE
  noise_init(lat, -1);
#endif
  
  return lat;
}

/************************************************************/

void d3q15_destroy(Lattice *lat) {
#ifdef DQ_NOISE
  noise_del(lat);
#endif

/*   free(lat->w); */
/*   free(lat->xi); */

  free(lat->f_ptr);
  free(lat->rho_ptr);
  free(lat->u_ptr);
  free(lat->force_ptr);
  
  free(lat);
}

/************************************************************/

void d3q15_iterate(Lattice *lat, int n_steps) {
  
  
  for (int t=0; t<n_steps; t++) {
    /* We have to work out the force first */
    (*lat->force_func)(lat);
    
    /* Collide */
    collide(lat);
    
    /* Update boundary */
    (*lat->bc_func)(lat);
    
    /* Propagate */
    propagate(lat);
    
    lat->time_step++;
  }
}

/************************************************************/

void propagate (Lattice *lat) {
  /* The propagation step of the algorithm.
   *
   * Copies (f[x][y][z][i] + R[i]) into f[x+dx][y+dy][z+dz][i]
   *
   * For each site, only want to access cells which are further
   * along in memory, so use: (expecting 7 directions, as
   * zero velocity doesn't move and we're swapping)
   * [1,0,0] [0,1,0] [0,0,1] i.e. 1, 3, 5
   * [1,1,1] [1,1,-1] [1,-1,1] [1,-1,-1] i.e. 7,8,9,10
   *
   * with, respectively:
   * 2,4,6
   * 14, 13, 12, 11
   *
   * We swap the value with the opposite direction velocity in 
   * the target Lattice site. By starting in the halo, we ensure
   * all the real cells get fully updated. Note that the first 
   * row must be treated a little differently, as it has no
   * neighbour in the [1,-1] position.
   * After this is complete, we have the distribution functions
   * updated, but with the fluid propagating AGAINST the velocity
   * vector. So need to reorder the values once we've done the 
   * propagation step.
   */
  double tmp;
  int i, j, k;
  int nx = lat->nx, ny = lat->ny, nz = lat->nz;
  /* Special cases POSSIBLY needed are the eight corners,
   * twelve edges and six faces. */


  for (i=0; i<=nx+1; i++) {
    for (j=0; j<=ny+1; j++) {
      for (k=0; k<=nz+1; k++) {
	/* propagate */
	
	/* [1,0,0] */
	if (i<=nx)
	  swap(DQ_f_get(lat, i,j,k, 1), DQ_f_get(lat, i+1,j,k, 2), tmp);
	/* [0,1,0] */
	if (j<=ny)
	  swap(DQ_f_get(lat, i,j,k, 3), DQ_f_get(lat, i,j+1,k, 4), tmp);
	/* [0,0,1] */
	if (k<=nz)
	  swap(DQ_f_get(lat, i,j,k, 5), DQ_f_get(lat, i,j,k+1, 6), tmp);
	
	/* [1,1,1] */
	if (i<=nx && j<=ny && k<= nz) 
	  swap(DQ_f_get(lat, i,j,k, 7), DQ_f_get(lat, i+1, j+1, k+1, 14), tmp);
	/* [1,1,-1] */
	if (i<=nx && j<=ny && k>0)
	  swap(DQ_f_get(lat, i,j,k, 8), DQ_f_get(lat, i+1, j+1, k-1, 13), tmp);
	/* [1,-1,1] */
	if (i<=nx && j>0 && k<=nz)
	  swap(DQ_f_get(lat, i,j,k, 9), DQ_f_get(lat, i+1, j-1, k+1, 12), tmp);
	/* [1,-1,-1] */
	if (i<=nx && j>0 && k >0)
	  swap(DQ_f_get(lat, i,j,k, 10),DQ_f_get(lat, i+1, j-1, k-1, 11), tmp);
	
	/* reorder */
	//if (i<=nx)
	swap(DQ_f_get(lat, i,j,k, 1), DQ_f_get(lat, i,j,k, 2), tmp);
	//if (j<=ny)
	swap(DQ_f_get(lat, i,j,k, 3), DQ_f_get(lat, i,j,k, 4), tmp);
	//if (k<=nz)
	swap(DQ_f_get(lat, i,j,k, 5), DQ_f_get(lat, i,j,k, 6), tmp);
	
	//if (i<=nz && j<=ny && k<=nz)
	swap(DQ_f_get(lat, i,j,k, 7), DQ_f_get(lat, i,j,k, 14), tmp);
	//if (i<=nz && j<=ny && k>0)
	swap(DQ_f_get(lat, i,j,k, 8), DQ_f_get(lat, i,j,k, 13), tmp);
	//if (i<=nz && j>0 && k<=nz)
	swap(DQ_f_get(lat, i,j,k, 9), DQ_f_get(lat, i,j,k, 12), tmp);
	//if (i<=nz && j>0 && k>0)
	swap(DQ_f_get(lat, i,j,k, 10),DQ_f_get(lat, i,j,k, 11), tmp);
      }
    }
  }
}

/************************************************************/

void calc_hydro_site(Site *site, Lattice* lat) {
  int i, a;
  double rho = 0;
  double mom[DQ_d];
  for (a=0; a<DQ_d; a++) {
    mom[a] = site->force[a] / 2.0;
  }
  
  for (i=0; i<DQ_q; i++) {
    rho += site->f[i];
    for (a=0; a<DQ_d; a++) {
      mom[a] += site->f[i];
    }
  }
  site->rho[0] = rho;
  for (a=0; a<DQ_d; a++) {
    site->u[a] = mom[a] / rho;
  }
}

void collide (Lattice *lat) {
  /* loop over cells */
  int i,j,k;
  /* loop indices for dimension */
  int a,b;
  /* loop indices for velocities & modes */
  int p;
  
  double fEq[DQ_q];
  double phi[DQ_q];
  Site site;
  
  double tau_s = lat->tau_s;
  double omega_s = 1.0 / (tau_s + 0.5);
  
  for (i=1; i<=lat->nx; i++) {
    for (j=1; j<=lat->ny; j++) {
      for (k=1; k<=lat->nz; k++) {
	set_site(lat, site, i,j,k);

	calc_hydro_site(&site, lat);
	// rho & u evaluated at t
	calc_equil(lat, site.rho[0], site.u, fEq);

	for (p=0; p<DQ_q; p++) {
	  // forcing (Phi) term PRE Eq. (5)
	  phi[p] = 0.0;
	  for (a=0; a<DQ_d; a++) {
	    phi[p] += site.force[a] * lat->xi[p][a] / lat->cs2;
	    for(b=0; b<DQ_d; b++) {
	      phi[p] += (site.u[a]*site.force[b] + site.u[b]*site.force[a])* lat->Q[p][a][b]
		/ (2.0 * lat->cs2 * lat->cs2 * site.rho[0]);
	    }
	  }
	  phi[p] *= lat->w[p];

	  /* Collide - Eq. (17) */
	  site.f[p] += omega_s * (tau_s * phi[p] - (site.f[p] - fEq[p]));
	}
      }	/* k */
    } /* j */
  } /* i */
  
}

/************************************************************/

void calc_hydro(Lattice *lat) {
  /* Work out the hydrodynamic variables for the whole grid. */
  int i,j,k;
  Site site;
  
  for (i=1; i<=lat->nx; i++) {
    for (j=1; j<=lat->ny; j++) {
      for (k=1; k<=lat->nz; k++) {
	set_site(lat, site, i,j,k);
	calc_hydro_site(&site, lat);
      }
    }
  }
}

/************************************************************/

void calc_equil(Lattice *lat, double rho, double u[], double f_eq[]) {
  int i, a, b;
  const double cs2 = lat->cs2;
  
  for (i = 0; i< DQ_q; i++) {
    double feqi = 1.0;
    for (a = 0; a < DQ_d; a++) {
      feqi += u[a]*lat->xi[i][a] / cs2;

      for (b = 0; b < DQ_d; b++) {
	feqi += lat->Q[i][a][b] * u[a] * u[b] / (2 * cs2 * cs2);
      }
    }
    f_eq[i] = rho * lat->w[i] * feqi;
  }
}

/************************************************************/

void calc_phi(double force[], double rho, double u[], double ans[]) {
  /* Calculate Phi, given by:
   * Phi_i = -rho*w_i * (F.xi_i/cs2 + (uF + Fu)Q_i/2cs4)
   *       =  rho*w_i * (F.u/cs2 - F.xi_i/cs2*(1+u.xi_i/cs2))
   */
  double rho_w0 = rho * 2.0 / 9.0,
    rho_w1 = rho / 9.0,
    rho_w2 = rho / 72.0;
  double cs2 = 1.0/ 3.0;
  double Fx_cs2 = force[0]/cs2, Fy_cs2 = force[1]/cs2, Fz_cs2 = force[2]/cs2;
  double Fxu_cs2 = Fx_cs2*u[0], Fyv_cs2 = Fy_cs2*u[1], Fzw_cs2 = Fz_cs2*u[2];
  
  double Fdotu_cs2 = Fxu_cs2 + Fyv_cs2 + Fzw_cs2;
  double u_cs2 = u[0]/cs2,      v_cs2 = u[1]/cs2,      w_cs2 = u[2]/cs2;

  ans[0] = rho_w0 * (Fdotu_cs2                                              );
  
  ans[1] = rho_w1*(Fdotu_cs2-( Fx_cs2              )*(1.0+u_cs2            ));
  ans[2] = rho_w1*(Fdotu_cs2-(-Fx_cs2              )*(1.0-u_cs2            ));
  ans[3] = rho_w1*(Fdotu_cs2-(        Fy_cs2       )*(1.0      +v_cs2      ));
  ans[4] = rho_w1*(Fdotu_cs2-(       -Fy_cs2       )*(1.0      -v_cs2      ));
  ans[5] = rho_w1*(Fdotu_cs2-(               Fz_cs2)*(1.0            +w_cs2));
  ans[6] = rho_w1*(Fdotu_cs2-(              -Fz_cs2)*(1.0            -w_cs2));
  
  ans[ 7] = rho_w2*(Fdotu_cs2-(+Fx_cs2+Fy_cs2+Fz_cs2)*(1.0+u_cs2+v_cs2+w_cs2));
  ans[ 8] = rho_w2*(Fdotu_cs2-(+Fx_cs2+Fy_cs2-Fz_cs2)*(1.0+u_cs2+v_cs2-w_cs2));
  ans[ 9] = rho_w2*(Fdotu_cs2-(+Fx_cs2-Fy_cs2+Fz_cs2)*(1.0+u_cs2-v_cs2+w_cs2));
  ans[10] = rho_w2*(Fdotu_cs2-(+Fx_cs2-Fy_cs2-Fz_cs2)*(1.0+u_cs2-v_cs2-w_cs2));
  ans[11] = rho_w2*(Fdotu_cs2-(-Fx_cs2+Fy_cs2+Fz_cs2)*(1.0-u_cs2+v_cs2+w_cs2));
  ans[12] = rho_w2*(Fdotu_cs2-(-Fx_cs2+Fy_cs2-Fz_cs2)*(1.0-u_cs2+v_cs2-w_cs2));
  ans[13] = rho_w2*(Fdotu_cs2-(-Fx_cs2-Fy_cs2+Fz_cs2)*(1.0-u_cs2-v_cs2+w_cs2));
  ans[14] = rho_w2*(Fdotu_cs2-(-Fx_cs2-Fy_cs2-Fz_cs2)*(1.0-u_cs2-v_cs2-w_cs2));
  
}

/************************************************************/

