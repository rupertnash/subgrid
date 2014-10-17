#ifndef DQ_PRNG
#define DQ_PRNG

#define DQ_PRNG_NTAB 32
typedef struct mt_state {
  long idum;
  long iy;
  long iv[DQ_PRNG_NTAB];
} mt_state;

void mt_init(mt_state* st, long idum);
void mt_del(mt_state* st);
double mt_get(mt_state* st);

typedef struct gasdev_state {
  int iset;
  double gset;
  mt_state mt_st;
} gasdev_state;

void gasdev_init(gasdev_state* gds, long idum);
void gasdev_del(gasdev_state* gds);
double gasdev_get(gasdev_state* gds);

typedef struct pgasdev_state {
  int n;
  gasdev_state* ptrs;
} pgasdev_state;

void pgasdev_init(pgasdev_state* pstate, int nThreads, long idum);
void pgasdev_del(pgasdev_state* gds);
double pgasdev_get(pgasdev_state* gds);

#endif
