#include "prng.h"
#include "omp.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void pgasdev_init(pgasdev_state* pgds, int nThreads, long idum) {
  if (nThreads <= 0)
    nThreads = omp_get_max_threads();
  
  pgds->n = nThreads;
  pgds->ptrs = (gasdev_state*)malloc(nThreads*sizeof(gasdev_state));
  
  for (int i=0; i<nThreads; i++) {
    gasdev_init(pgds->ptrs + i, idum - i);
  }
}
void pgasdev_del(pgasdev_state* pgds) {
  for (int i=0; i<pgds->n; i++) {
    gasdev_del(pgds->ptrs + i);
  }
  free(pgds->ptrs);
}

double pgasdev_get(pgasdev_state* pgds) {
  int i = omp_get_thread_num();
  return gasdev_get(pgds->ptrs + i);
}

void gasdev_init(gasdev_state* gds, long idum) {
  gds->iset = 0;
  mt_init(&gds->mt_st, idum);
}

void gasdev_del(gasdev_state* gds) {
  mt_del(&gds->mt_st);
}

double gasdev_get(gasdev_state* gds)
     /* from Numerical recipes */
{
    double fac,rsq,v1,v2;

    if (gds->iset == 0) {
        do {
            v1=2.0*mt_get(&gds->mt_st)-1.0;
            v2=2.0*mt_get(&gds->mt_st)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gds->gset=v1*fac;
        gds->iset=1;
        return (double)(v2*fac);
    } else {
        gds->iset=0;
        return (double)gds->gset;
    }
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/DQ_PRNG_NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


void mt_init(mt_state* st, long idum) {
  int j;
  long k;
  if (-idum < 1)
    idum=1;
  else
    idum = -idum;

  for (j=DQ_PRNG_NTAB+7; j>=0; j--) {
    k= idum/IQ;
    idum = IA*(idum - k*IQ)-IR*k;
    if (idum < 0)
      idum += IM;
    if (j < DQ_PRNG_NTAB)
      st->iv[j] = idum;
  }
  st->iy=st->iv[0];
  st->idum = idum;
}

void mt_del(mt_state* st) {
}

double mt_get(mt_state* st) {
    int j;
    long k;
    double temp;
    
    if (st->idum <= 0 || !st->iy) {
      printf("Twister error\n");
      exit(1);
    }
    k=(st->idum)/IQ;
    st->idum=IA*(st->idum-k*IQ)-IR*k;
    if (st->idum < 0) st->idum += IM;
    j=st->iy/NDIV;
    st->iy=st->iv[j];
    st->iv[j] = st->idum;
    temp=AM*st->iy;
    if (temp > RNMX) return RNMX;
    else return temp;
}
