#ifndef D3Q15_H
#define D3Q15_H

#define DQ_NEW 1
/* Dimension */
#define DQ_d 3
/* number of velocities */
#define DQ_q 15


/* Macro to swap values */
#define swap(a,b,temp) {temp = a; a = b; b = temp;}



/* Macros for accessing  elements in the way you'd expect. */
#define DQ_f_get(L, i,j,k, m) *(L->f_current_ptr + ((i) * L->strides[DQ_X] +\
					 (j) * L->strides[DQ_Y] +\
					 (k) * L->strides[DQ_Z]) * DQ_q + (m))
#define DQ_f_new_get(L, i,j,k, m) *(L->f_new_ptr + ((i) * L->strides[DQ_X] +\
					 (j) * L->strides[DQ_Y] +\
					 (k) * L->strides[DQ_Z]) * DQ_q + (m))

/* These could easily be replaced with functions to save memory */
#define DQ_rho_get(L, i,j,k) *(L->rho_ptr + ((i) * L->strides[DQ_X] +\
					  (j) * L->strides[DQ_Y] + \
					  (k) * L->strides[DQ_Z]))
#define DQ_u_get(L, i,j,k, m) *(L->u_ptr + ((i) * L->strides[DQ_X] + \
					 (j) * L->strides[DQ_Y] + \
					 (k) * L->strides[DQ_Z]) * DQ_d + (m))
#define DQ_force_get(L, i,j,k, m) *(L->force_ptr + ((i) * L->strides[DQ_X] + \
						 (j) * L->strides[DQ_Y] + \
						 (k) * L->strides[DQ_Z])*DQ_d+(m))
#define DQ_access_macros
#include "eigenvectors.h"
#undef DQ_access_macros

typedef struct NoiseConfig NoiseConfig;
#include "Lattice.h"
#include "noise.h"

/* function declarations */
Lattice *d3q15_init(int nx, int ny, int nz, double tau_s, double tau_b);
void d3q15_iterate(Lattice *lat, int n_steps);
void d3q15_step(Lattice *lat);
void d3q15_destroy(Lattice *lat);

void dq_push(Lattice *lat, const int i, const int j, const int k, const double fPostCollision[DQ_q]);
void dq_collide(const Lattice *lat, const int i, const int j, const int k, double fPostCollision[DQ_q]);

void calc_equil(double rho, double *u, double *f_eq);
void calc_hydro(Lattice *lat);

/* Forces... */
void force_none_init(Lattice *lat);
void force_none_calc(Lattice *lat);

void force_const_init(Lattice *lat, double F[]);
void force_const_calc(Lattice *lat);

/* Boundary conditions */
void bc_pbc_init(Lattice *lat);
void bc_pbc_update(Lattice *lat);
void bc_pbc_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d]);

void bc_noslip_init(Lattice *lat);
void bc_noslip_update(Lattice *lat);
void bc_noslip_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d]);

void bc_freeslip_init(Lattice *lat);
void bc_freeslip_update(Lattice *lat);
void bc_freeslip_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d]);

void bc_wall_init(Lattice *lat,
		  double v_up_amp[DQ_d], double v_up_off[DQ_d],
		  double v_lo_amp[DQ_d], double v_lo_off[DQ_d],
		  double shear_freq, double shear_phase);
void bc_wall_update(Lattice *lat);
void bc_wall_do_site(Lattice *lat, int i[DQ_d], int n[DQ_d], double v[DQ_d]);

void bc_pressure_init(Lattice *lat);
void bc_pressure_update(Lattice *lat);

/* Initial conditions */
void init_uniform (Lattice *lat, double rho, double u[]);

void init_point (Lattice *lat, double rho, double u[],
		 double rho_prime, int x[]);

/* Checks */
void total_mass_and_momentum(Lattice *lat, double *mass, double mom[DQ_d]);

#endif
