#ifndef DQ_NOISE
#define DQ_NOISE

struct NoiseConfig {
  double var[5];
  double _var[5];
  double temperature;
  double rootT;
  
  long seed;
};

void noise_init(Lattice *lat, long seed);
void noise_set_temperature(NoiseConfig *noise, double temperature);
void noise_add_to_modes(NoiseConfig *n, double mode[]);
void noise_del(Lattice *lat);

double gasdev(long *idum);
double ran1(long *idum);

#endif
