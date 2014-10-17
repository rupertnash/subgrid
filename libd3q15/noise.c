#include "d3q15.h"

#ifdef DQ_NOISE

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void noise_init(Lattice *lat, long seed) {
  NoiseConfig *n = (NoiseConfig *)malloc(sizeof(NoiseConfig));

  double omega_s = 1.0 / (lat->tau_s + 0.5);
  double omega_b = 1.0 / (lat->tau_b + 0.5);
  
  pgasdev_init(&n->gds, -1, seed);
  
  n->_var[0] = 
    sqrt(1./lat->norms[DQ_SXX]) *sqrt(2.* lat->tau_b * omega_b*omega_b);
  n->_var[1] = 
    sqrt(1./lat->norms[DQ_SXY]) *sqrt(2.* lat->tau_s * omega_s*omega_s);
  
  /* Multirelaxtion sets tau = 1 in the final three variances.
   * Note that Ladd's has these three zero. */
  n->_var[2] = sqrt(1./lat->norms[DQ_chi1]);
  n->_var[3] = sqrt(1./lat->norms[DQ_jchi1X]);
  n->_var[4] = sqrt(1./lat->norms[DQ_chi2]);
  
  noise_set_temperature(n, 0.0);
  
  lat->noise = n;
}

/************************************************************/

void noise_set_temperature(NoiseConfig *noise, double temperature){
  int i;
  noise->temperature = temperature;
  noise->rootT = sqrt(temperature);
  for (i=0; i<5; i++)
    noise->var[i] = noise->rootT*noise->_var[i];
  
}

/************************************************************/

void noise_del(Lattice *lat) {
  pgasdev_del(&lat->noise->gds);
  free(lat->noise);
}

/************************************************************/

void noise_add_to_modes(NoiseConfig *n, double mode[]) {
  int a;
  double Shat[DQ_d][DQ_d];
  double TrS;
  double chi1hat;
  double jchi1hat[DQ_d];
  double chi2hat;
  
  /* density & momentum modes unchanged */

  /* stress mode */
  Shat[DQ_X][DQ_X] = pgasdev_get(&n->gds);
  Shat[DQ_X][DQ_Y] = pgasdev_get(&n->gds);
  Shat[DQ_X][DQ_Z] = pgasdev_get(&n->gds);

  Shat[DQ_Y][DQ_X] = Shat[DQ_X][DQ_Y];
  Shat[DQ_Y][DQ_Y] = pgasdev_get(&n->gds);
  Shat[DQ_Y][DQ_Z] = pgasdev_get(&n->gds);

  Shat[DQ_Z][DQ_X] = Shat[DQ_X][DQ_Z];
  Shat[DQ_Z][DQ_Y] = Shat[DQ_Y][DQ_Z];
  Shat[DQ_Z][DQ_Z] = pgasdev_get(&n->gds);
  
  /* compute trace & traceless */
  TrS = 0.;
  for (a=0; a<DQ_d; a++)
    TrS += Shat[a][a];
  TrS /= 3.;
  
  for (a=0; a<DQ_d; a++)
    Shat[a][a] -= TrS;
  
  /* set variances */
  Shat[DQ_X][DQ_X] *= sqrt(2.) * n->var[1];
  Shat[DQ_X][DQ_Y] *= n->var[1];
  Shat[DQ_X][DQ_Z] *= n->var[1];
  
  Shat[DQ_Y][DQ_X] *= n->var[1];
  Shat[DQ_Y][DQ_Y] *= sqrt(2.) * n->var[1];
  Shat[DQ_Y][DQ_Z] *= n->var[1];

  Shat[DQ_Z][DQ_X] *= n->var[1];
  Shat[DQ_Z][DQ_Y] *= n->var[1];
  Shat[DQ_Z][DQ_Z] *= sqrt(2.) * n->var[1];
  
  TrS *= n->var[0];
  
  /* add trace back on */
  for (a=0; a<DQ_d; a++)
    Shat[a][a] += TrS;
  
  chi1hat      = n->var[2] * pgasdev_get(&n->gds);
  
  jchi1hat[0]  = n->var[3] * pgasdev_get(&n->gds);
  jchi1hat[1]  = n->var[3] * pgasdev_get(&n->gds);
  jchi1hat[2]  = n->var[3] * pgasdev_get(&n->gds);
  
  chi2hat      = n->var[4] * pgasdev_get(&n->gds);

  mode[DQ_SXX] += Shat[DQ_X][DQ_X];
  mode[DQ_SXY] += Shat[DQ_X][DQ_Y];
  mode[DQ_SXZ] += Shat[DQ_X][DQ_Z];
  mode[DQ_SYY] += Shat[DQ_Y][DQ_Y];
  mode[DQ_SYZ] += Shat[DQ_Y][DQ_Z];
  mode[DQ_SZZ] += Shat[DQ_Z][DQ_Z];
  
  mode[DQ_chi1] += chi1hat;
  mode[DQ_jchi1X] += jchi1hat[DQ_X];
  mode[DQ_jchi1Y] += jchi1hat[DQ_Y];
  mode[DQ_jchi1Z] += jchi1hat[DQ_Z];

  mode[DQ_chi2] += chi2hat;

}

/************************************************************/

void noise_calc(Lattice *lat, double fr[]) {
  int i, j, p;
  
  double sigma[3][3];       /* Random (symmetric) matrix */
  double jchihat[3];        /* Random vector */
  double chi1hat;           /* Random */
  double chi2hat;           /* Random */

  double var;
  double rhat[11];

  NoiseConfig *n = lat->noise;

  double *_var = n->var; /* pointer to variances */

  for (i = 0; i < 11; i++) {
    rhat[i] = pgasdev_get(&n->gds);
  }

  sigma[0][0] = _var[0]*rhat[0];
  sigma[1][1] = _var[0]*rhat[1];
  sigma[2][2] = _var[0]*rhat[2];

  sigma[0][1] = _var[1]*rhat[3];
  sigma[1][2] = _var[1]*rhat[4];
  sigma[2][0] = _var[1]*rhat[5];

  sigma[1][0] = sigma[0][1];
  sigma[2][1] = sigma[1][2];
  sigma[0][2] = sigma[2][0];

  chi1hat     = _var[2]*rhat[6];

  jchihat[0]  = _var[3]*rhat[7];
  jchihat[1]  = _var[3]*rhat[8];
  jchihat[2]  = _var[3]*rhat[9];

  chi2hat     = _var[4]*rhat[10];

  for (p = 0; p < DQ_q; p++) {
    var = 0.0;
    for (i = 0; i < DQ_d; i++) {
      for (j = 0; j < DQ_d; j++) {
	var += sigma[i][j]*lat->Q[p][i][j];
      }
    }

    fr[p]  = 4.5*var;
    fr[p] += 0.5*chi1hat* lat->mm[DQ_chi1][p];

    var = 0.0;
    for (i = 0; i < DQ_d; i++) {
      var += jchihat[i]* lat->mm[DQ_jchi1X + i][p];
    }
    
    fr[p] += 1.5*var;
    fr[p] += 9. * lat->mm[DQ_chi2][p] *chi2hat;
    fr[p] *= lat->w[p];
  }

}

/************************************************************/



#endif
