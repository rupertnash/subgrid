#include "CDelta.h"
#include <math.h>

double CDelta_delta(double x) {
  double abs_x = fabs(x);
  double root = -4. * x*x;
  double phi = -2.* abs_x;
  
  if (abs_x >= 2.0)
    return 0.;
  
  if (abs_x >= 1.0) {
    root += 12. * abs_x - 7.;
    phi += 5.;
    phi -= sqrt(root);
  } else {
    root += 4. * abs_x + 1;
    phi += 3.;
    phi += sqrt(root);
  }
  return 0.125 * phi;
}

/* double CDelta_delta(double x) { */
/*   double abs_x = fabs(x); */
/*   if (abs_x >= 2.0) { */
/*     return 0.; */
/*   } else if (abs_x >= 1.0) { */
/*     return (5. - 2.*abs_x - sqrt(-7. + 12.*abs_x - 4.*abs_x*abs_x) ) * 0.125; */
/*   } else { */
/*     return (3. - 2.*abs_x + sqrt( 1. +  4.*abs_x - 4.*abs_x*abs_x) ) * 0.125; */
/*   } */
/* } */

double CDelta_3d_delta(double x[3]) {
  double ans = 1.;
  int d;
  for (d=0; d<3; d++) {
    ans *= CDelta_delta(x[d]);
  }
  return ans;
}
