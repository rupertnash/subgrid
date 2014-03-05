#ifndef DQ_Lattice_H
#define DQ_Lattice_H
/* A Lattice type for ease of passing the important parameters 
 * and (pointers to) the Lattice around.
 */
typedef struct Lattice {
  /* Could be useful anywhere */

  /* speed of sound, squared */
  double cs2;
  /* quadrature weights */
  double w[DQ_q];
  /* velocities */
  double xi[DQ_q][DQ_d];
  /* Kinetic projector, Q_i_ab = xi_i_a * xi_i_b - delta_ab cs^2 */
  double Q[DQ_q][DQ_d][DQ_d];
  
  int complement[DQ_q];
  
  /* The modes' normalizers */
  double norms[DQ_q];
  /* Matrix of eigenvectors (they are the rows) */
  double mm[DQ_q][DQ_q];
  /* and its inverse */
  double mmi[DQ_q][DQ_q];
  
  /* Size of the Lattice to use */
  int nx;
  int ny;
  int nz;

  int strides[DQ_d];
  
  /* Size of the arrays used to store Lattice,
   * needed for visualisation with libgraph */
  int _x_ar_size;
  int _y_ar_size;
  int _z_ar_size;
  
  /* Pointers to the distribution function array and 
   * hydrodynamic variables (so that they can be output/visualised 
   * more easily. 
   */
  double *f_ptr;
  double *rho_ptr;
  double *u_ptr;
  double *force_ptr;
  
  /* The relaxation time */
  double tau_s; /* shear relaxation time */
  double tau_b; /* bulk relaxation time */
  
  /* Current timestep */
  int time_step;
  
  /* Pointers to functions which do configuation dependent stuff */
  void (*bc_func)(struct Lattice *);
  void (*force_func)(struct Lattice *);

  /* Pointer to the noise config. */
  NoiseConfig *noise;

  /* A pointer to whatever other information might be needed by 
   * the controlling code. Casts to and from void will be needed.
   */
  void *config;
} Lattice;

typedef struct Site {
  double *f;
  double *rho;
  double *u;
  double *force;
} Site;

#define set_site(L, s, i,j,k) {int ijk= ((i) * L->strides[DQ_X] + \
			  	         (j) * L->strides[DQ_Y] + \
				         (k) * L->strides[DQ_Z]); \
                               s.f     = L->f_ptr     + ijk*DQ_q; \
                               s.rho   = L->rho_ptr   + ijk;	  \
                               s.u     = L->u_ptr     + ijk*DQ_d; \
                               s.force = L->force_ptr + ijk*DQ_d; }

#endif
