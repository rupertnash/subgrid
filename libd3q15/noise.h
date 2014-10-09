#ifndef DQ_NOISE
#define DQ_NOISE

typedef struct gasdev_state gasdev_state;

struct NoiseConfig {
  double var[5];
  double _var[5];
  double temperature;
  double rootT;
  
  gasdev_state* gds;
};

void noise_init(Lattice *lat, long seed);
void noise_set_temperature(NoiseConfig *noise, double temperature);
void noise_add_to_modes(NoiseConfig *n, double mode[]);
void noise_del(Lattice *lat);

gasdev_state* gasdev_init(long idum);
double gasdev_get(gasdev_state* gds);
void gasdev_del(gasdev_state* gds);

#endif
