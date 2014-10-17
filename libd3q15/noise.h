#ifndef DQ_NOISE
#define DQ_NOISE

#include "prng.h"

struct NoiseConfig {
  double var[5];
  double _var[5];
  double temperature;
  double rootT;
  
  pgasdev_state gds;
};

void noise_init(Lattice *lat, long seed);
void noise_set_temperature(NoiseConfig *noise, double temperature);
void noise_add_to_modes(NoiseConfig *n, double mode[]);
void noise_del(Lattice *lat);


#endif
