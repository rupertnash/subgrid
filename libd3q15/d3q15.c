/* 3 dimensional, 15 velocity Lattice Boltzmann code.
 */
#include <stdlib.h>
#include <math.h>
#include "d3q15.h"
#include "omp.h"

/************************************************************/

/* static const int q = DQ_q; */
/* static const int d = DQ_d; */

/* identity matrix */
static double DQ_delta[DQ_d][DQ_d];
static int DQ_delta_initialised = 0;

#define DQ_MAX_BLOCK_SIZE 8
typedef struct Block {
  int min[DQ_d];
  int max[DQ_d];
} Block;
struct BlockSched {
  int n;
  Block* blocks;
};


Lattice *d3q15_init(int nx, int ny, int nz, double tau_s, double tau_b) {
  if (!DQ_delta_initialised) {
    for (int a=0; a<DQ_d; a++) {
      for (int b=0; b<DQ_d; b++) {
	DQ_delta[a][b] = 0.;
      }
      DQ_delta[a][a] = 1./DQ_d;
    }
    DQ_delta_initialised = 1;
  }
  
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
  lat->f_current_ptr = (double *)malloc( (nx+2) * (ny+2) * (nz+2) * DQ_q *
					 sizeof(double) );
  lat->f_new_ptr = (double *)malloc( (nx+2) * (ny+2) * (nz+2) * DQ_q *
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
  
  int nSites[DQ_d] = {nx, ny, nz};
  int nBlocks[DQ_d];
  int totalBlocks = 1;
  int* bounds[DQ_d];
  for (int d=0; d<DQ_d; d++) {
    nBlocks[d] = (int)ceil((double)nSites[d] / DQ_MAX_BLOCK_SIZE);
    totalBlocks *= nBlocks[d];
    double delta = (double)nSites[d] / nBlocks[d];
    bounds[d] = (int*)malloc((nBlocks[d]+1)*sizeof(int));
    for (int i=0; i<=nBlocks[d]; i++)
      bounds[d][i] = (int)(delta * i) + 1;
  }
  
  lat->block_schedule = (BlockSched*)malloc(sizeof(BlockSched));
  lat->block_schedule->n = totalBlocks;
  lat->block_schedule->blocks = (Block*)malloc(totalBlocks*sizeof(Block));
  
  int iBlock[DQ_d];
  int i = 0;
  for (iBlock[0]=0; iBlock[0]<nBlocks[0]; iBlock[0]++) {
  for (iBlock[1]=0; iBlock[1]<nBlocks[1]; iBlock[1]++) {
  for (iBlock[2]=0; iBlock[2]<nBlocks[2]; iBlock[2]++) {
    for (int d=0; d<DQ_d; d++) {
      lat->block_schedule->blocks[i].min[d] = bounds[d][iBlock[d]];
      lat->block_schedule->blocks[i].max[d] = bounds[d][iBlock[d]+1];
    }
    i++;
  }}}
  for (int d=0; d<DQ_d; d++)
    free(bounds[d]);

  /* Now, go and touch the f's from the core that will be in charge of
     them to ensure they're allocated in the proper NUMA region. */
#pragma omp parallel for
  for (int iBlock=0; iBlock<lat->block_schedule->n; iBlock++) {
    
    int min[DQ_d];
    int max[DQ_d];
    for (int d=0; d<DQ_d; d++) {
      min[d] = lat->block_schedule->blocks[iBlock].min[d];
      min[d] = (min[d] == 1 ? 0 : min[d]);
      max[d] = lat->block_schedule->blocks[iBlock].max[d];
      max[d] = (max[d] == nSites[d] ? nSites[d]+1 : max[d]);
    }
    for (int i=min[0]; i<max[0]; i++) {
    for (int j=min[1]; j<max[1]; j++) {
    for (int k=min[2]; k<max[2]; k++) {
      for (int p=0; p<DQ_q; p++) {
	DQ_f_get(lat, i, j, k, p) = 0.0;
	DQ_f_new_get(lat, i, j, k, p) = 0.0;
      }
    }}}
  }
  return lat;
}

/************************************************************/

void d3q15_destroy(Lattice *lat) {
#ifdef DQ_NOISE
  noise_del(lat);
#endif
  free(lat->block_schedule);
  free(lat->f_current_ptr);
  free(lat->f_new_ptr);
  free(lat->rho_ptr);
  free(lat->u_ptr);
  free(lat->force_ptr);
  
  free(lat);
}

/************************************************************/

void d3q15_iterate(Lattice *lat, int n_steps) {
  for (int t=0; t<n_steps; t++) {
    d3q15_step(lat);
  }
}

void d3q15_step(Lattice* lat) {
  /* We have to work out the force first */
  (*lat->force_func)(lat);

#pragma omp parallel for
  for (int iBlock=0; iBlock<lat->block_schedule->n; iBlock++) {
    
    int* min = lat->block_schedule->blocks[iBlock].min;
    int* max = lat->block_schedule->blocks[iBlock].max;
    
    for (int i=min[0]; i<max[0]; i++) {
    for (int j=min[1]; j<max[1]; j++) {
    for (int k=min[2]; k<max[2]; k++) {
      double fPostCollision[DQ_q];
      /* Collide */
      dq_collide(lat, i, j, k, fPostCollision);
      /* Propagate */
      dq_push(lat, i, j ,k, fPostCollision);
    }}}
  }
  
  /* Update boundary */
  (*lat->bc_func)(lat);
  
  /* Swap the distribution arrays */
  double* tmp = lat->f_current_ptr;
  lat->f_current_ptr = lat->f_new_ptr;
  lat->f_new_ptr = tmp;
  
  lat->time_step++;
}

/************************************************************/

void dq_push(Lattice *lat, const int i, const int j, const int k, const double fPostCollision[DQ_q]) {
  /* Propagate distributions using a push */
  const int* cp;
  for (int p=0; p<DQ_q; p++) {
    cp = lat->ci[p];
    DQ_f_new_get(lat,
		 i + cp[0],
		 j + cp[1],
		 k + cp[2],
		 p) = fPostCollision[p];
  }
}
void dq_pull(const Lattice *lat, const int i, const int j, const int k, double fPreCollision[DQ_q]) {
  /* Propagate distributions using a pull */
  const int* cp;
  for (int p=0; p<DQ_q; p++) {
    cp = lat->ci[p];
    fPreCollision[p] = DQ_f_get(lat,
	     i - cp[0],
	     j - cp[1],
	     k - cp[2],
	     p);
  }
}

/************************************************************/

void dq_collide (const Lattice *lat, const int i, const int j, const int k, double fPostCollision[DQ_q]) {
  /* loop indices for dimension */
  int a,b;
  /* loop indices for velocities & modes */
  int p,m;
  
  double mode[DQ_q];
  /* convenience vars */
  double S[DQ_d][DQ_d];
  double u[DQ_d];
  double usq, TrS, uDOTf;
  
  double tau_s = lat->tau_s;
  double tau_b = lat->tau_b;
  double omega_s = 1.0 / (tau_s + 0.5);
  double omega_b = 1.0 / (tau_b + 0.5);

  const double* f_t = &DQ_f_get(lat, i,j,k, 0);
  const double* force = &DQ_force_get(lat, i,j,k, 0);
  
  for (m=0; m<DQ_q; m++) {
    /* compute the modes */
    mode[m] = 0.;
    for (p=0; p<DQ_q; p++) {
      mode[m] += f_t[p] * lat->mm[m][p];
    }
  }
  
  /* This was a bad idea! EITHER rho is set correctly cos we've
   *  already calculated the hydrodynamic fields and we're performing a
   *  redundent store operation OR we're setting rho the value it
   *  should have at the start of the time step. (And similarly for u)
   */
  //site.rho[0] = mode[DQ_rho];
  
  /* Work out the site fluid velocity
   *   rho*u= (rho*u') + F*\Delta t /2
   * and the coefficient for the momentum modes.
   *    = (rho*u') * F* \Delta t
   * (and u squared)
   */
  usq = 0.;
  uDOTf = 0.;
  for (a=0; a<DQ_d; a++) {
    u[a] = (mode[DQ_mom(a)] + 0.5*force[a]) / mode[DQ_rho];
    mode[DQ_mom(a)] += force[a];
    usq += u[a]*u[a];
    uDOTf += u[a]*force[a];
  }
  
  /* For unequal relax trace & traceless part at different rates.
   * Equilibrium trace = rho*usq */
  
  /* First copy the stress to a convenience var */
  S[DQ_X][DQ_X] = mode[DQ_SXX];
  S[DQ_X][DQ_Y] = mode[DQ_SXY];
  S[DQ_X][DQ_Z] = mode[DQ_SXZ];
  
  S[DQ_Y][DQ_X] = mode[DQ_SXY];
  S[DQ_Y][DQ_Y] = mode[DQ_SYY];
  S[DQ_Y][DQ_Z] = mode[DQ_SYZ];
  
  S[DQ_Z][DQ_X] = mode[DQ_SXZ];
  S[DQ_Z][DQ_Y] = mode[DQ_SYZ];
  S[DQ_Z][DQ_Z] = mode[DQ_SZZ];
  
  /* Form the trace part */
  TrS = 0.;
  for (a=0; a<DQ_d; a++) {
    TrS += S[a][a];
  }
  /* And the traceless part */
  for (a=0; a<DQ_d; a++) {
    S[a][a] -= TrS/DQ_d;
  }
  
  /* relax the trace */
  TrS -= omega_b*(TrS - mode[DQ_rho]*usq);
  /* Add forcing part to trace */
  TrS += 2.*omega_b*tau_b * uDOTf;
  
  /* and the traceless part */
  for (a=0; a<DQ_d; a++) {
    for (b=0; b<DQ_d; b++) {
      S[a][b] -= omega_s*(S[a][b] - 
			  mode[DQ_rho]*(u[a]*u[b] -usq*DQ_delta[a][b]));
      
      /* including traceless force */
      S[a][b] += 2.*omega_s*tau_s * (u[a]*force[b] + force[a]*u[b] - 2. * uDOTf * DQ_delta[a][b]);
    }
    /* add the trace back on */
    S[a][a] += TrS / DQ_d;
  }
  
  /* copy S back into modes[] */
  mode[DQ_SXX] = S[DQ_X][DQ_X];
  mode[DQ_SXY] = S[DQ_X][DQ_Y];
  mode[DQ_SXZ] = S[DQ_X][DQ_Z];
	
  mode[DQ_SYY] = S[DQ_Y][DQ_Y];
  mode[DQ_SYZ] = S[DQ_Y][DQ_Z];
  
  mode[DQ_SZZ] = S[DQ_Z][DQ_Z];
  
  /* Ghosts are relaxed to zero immediately */
  mode[DQ_chi1] = 0.;
  mode[DQ_jchi1X] = 0.;
  mode[DQ_jchi1Y] = 0.;
  mode[DQ_jchi1Z] = 0.;
  mode[DQ_chi2] = 0.;
  
#ifdef DQ_NOISE
  noise_add_to_modes(lat->noise, mode);
#endif
  
  /* project back to the velocity basis */
  for (p=0; p<DQ_q; p++) {
   fPostCollision[p] = 0.;
    for (m=0; m<DQ_q; m++) {
      fPostCollision[p] += mode[m] * lat->mmi[p][m];
    }
  }
}

/************************************************************/

void calc_hydro(Lattice *lat) {
  /* Work out the hydrodynamic variables for the whole grid. */
  int i,j,k, m,p, a;
  double mode[DQ_q];
  Site site;
  
  for (i=1; i<=lat->nx; i++) {
    for (j=1; j<=lat->ny; j++) {
      for (k=1; k<=lat->nz; k++) {
	set_site(lat, site, i,j,k);
	
	for (m=0; m<DQ_q; m++) {
	  /* compute the modes */
	  mode[m] = 0.;
	  for (p=0; p<DQ_q; p++) {
	    mode[m] += site.f_current[p] * lat->mm[m][p];
	  }
	}
	
	site.rho[0] = mode[DQ_rho];
	
	/* Work out the site fluid velocity
	 *   rho*u= (rho*u') + F*\Delta t /2
	 */
	for (a=0; a<DQ_d; a++) {
	  site.u[a] = (mode[DQ_mom(a)] + 0.5*site.force[a]) / site.rho[0];
	}
	
      }
    }
  }
}

/************************************************************/

void calc_equil(double rho, double u[], double f_eq[]) {
  /* Works out the equilibrium distribution from the hydrodynamic
   * variables.
   */
  double cs2 = 1.0 / 3.0;
  double u2 = u[0]*u[0], v2 = u[1]*u[1], w2 = u[2]*u[2];
  double mod_sq = (u2 + v2 + w2)/(2.0 * cs2);
  double u_cs2 = u[0]/cs2, v_cs2 = u[1]/cs2, w_cs2 = u[2]/cs2;
  double u2_2cs4 = u2 / (2.0 * cs2 * cs2), v2_2cs4 = v2 / (2.0 * cs2 * cs2),
    w2_2cs4 = w2 / (2.0 * cs2 *cs2);
  double uv_cs4 = u_cs2*v_cs2, vw_cs4 = v_cs2*w_cs2, uw_cs4 = u_cs2*w_cs2;
  double mod_sq_2 = (u2 + v2 + w2) * (1 - cs2) / (2.0 * cs2 * cs2);

  double rho_w0 = rho * 2.0/9.0,
    rho_w1 = rho / 9.0,
    rho_w2 = rho / 72.0;
  
  f_eq[0] = rho_w0 * (1.0 - mod_sq);
  
  f_eq[1] = rho_w1 * (1.0 - mod_sq + u_cs2 + u2_2cs4);
  f_eq[2] = rho_w1 * (1.0 - mod_sq - u_cs2 + u2_2cs4);
  f_eq[3] = rho_w1 * (1.0 - mod_sq + v_cs2 + v2_2cs4);
  f_eq[4] = rho_w1 * (1.0 - mod_sq - v_cs2 + v2_2cs4);
  f_eq[5] = rho_w1 * (1.0 - mod_sq + w_cs2 + w2_2cs4);
  f_eq[6] = rho_w1 * (1.0 - mod_sq - w_cs2 + w2_2cs4);
  

  f_eq[7] = rho_w2 * (1.0 + u_cs2 + v_cs2 + w_cs2 +
		      uv_cs4 + vw_cs4 + uw_cs4 + mod_sq_2);
  
  f_eq[8] = rho_w2 * (1.0 + u_cs2 + v_cs2 - w_cs2 +
		      uv_cs4 - vw_cs4 - uw_cs4 + mod_sq_2);
  
  f_eq[9] = rho_w2 * (1.0 + u_cs2 - v_cs2 + w_cs2 -
		      uv_cs4 - vw_cs4 + uw_cs4 + mod_sq_2);
  
  f_eq[10] = rho_w2 * (1.0 + u_cs2 - v_cs2 - w_cs2 -
		      uv_cs4 + vw_cs4 - uw_cs4 + mod_sq_2);
  
  f_eq[11] = rho_w2 * (1.0 - u_cs2 + v_cs2 + w_cs2 -
		      uv_cs4 + vw_cs4 - uw_cs4 + mod_sq_2);
  
  f_eq[12] = rho_w2 * (1.0 - u_cs2 + v_cs2 - w_cs2 -
		      uv_cs4 - vw_cs4 + uw_cs4 + mod_sq_2);
  
  f_eq[13] = rho_w2 * (1.0 - u_cs2 - v_cs2 + w_cs2 +
		      uv_cs4 - vw_cs4 - uw_cs4 + mod_sq_2);
  
  f_eq[14] = rho_w2 * (1.0 - u_cs2 - v_cs2 - w_cs2 +
		      uv_cs4 + vw_cs4 + uw_cs4 + mod_sq_2);
  
}

/************************************************************/


