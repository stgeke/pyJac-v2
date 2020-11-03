#include "chem_utils.hpp"
#include <math.h>
#include <stdlib.h>




void eval_b(double const *__restrict__ phi, double *__restrict__ b)
{
  double T;
  double Tinv;
  double logT;

  {
    int const j = 0;

    logT = log(fmin(10000.0, fmax(100.0, phi[0])));
    Tinv = 1.0 / fmin(10000.0, fmax(100.0, phi[0]));
    T = fmin(10000.0, fmax(100.0, phi[0]));
    for (int k = 0; k <= 8; ++k)
    {
      if (T < T_mid[k])
        b[k] = T * (T * (T * (T * a_lo[4 + 7 * k] / 20.0 + a_lo[3 + 7 * k] / 12.0) + a_lo[2 + 7 * k] / 6.0) + a_lo[1 + 7 * k] / 2.0) + (a_lo[7 * k] + -1.0) * logT + -1.0 * a_lo[7 * k] + a_lo[6 + 7 * k] + -1.0 * a_lo[5 + 7 * k] * Tinv;
      if (!(T < T_mid[k]))
        b[k] = T * (T * (T * (T * a_hi[4 + 7 * k] / 20.0 + a_hi[3 + 7 * k] / 12.0) + a_hi[2 + 7 * k] / 6.0) + a_hi[1 + 7 * k] / 2.0) + (a_hi[7 * k] + -1.0) * logT + -1.0 * a_hi[7 * k] + a_hi[6 + 7 * k] + -1.0 * a_hi[5 + 7 * k] * Tinv;
    }
  }
}
void eval_h(double const *__restrict__ phi, double *__restrict__ h)
{
  double T;

  {
    int const j = 0;

    T = fmin(10000.0, fmax(100.0, phi[0]));
    for (int k = 0; k <= 8; ++k)
    {
      if (T < T_mid[k])
        h[k] = 8314.4621 * (T * (T * (T * (T * (T * a_lo[4 + 7 * k] / 5.0 + a_lo[3 + 7 * k] / 4.0) + a_lo[2 + 7 * k] / 3.0) + a_lo[1 + 7 * k] / 2.0) + a_lo[7 * k]) + a_lo[5 + 7 * k]);
      if (!(T < T_mid[k]))
        h[k] = 8314.4621 * (T * (T * (T * (T * (T * a_hi[4 + 7 * k] / 5.0 + a_hi[3 + 7 * k] / 4.0) + a_hi[2 + 7 * k] / 3.0) + a_hi[1 + 7 * k] / 2.0) + a_hi[7 * k]) + a_hi[5 + 7 * k]);
    }
  }
}
void eval_cp(double const *__restrict__ phi, double *__restrict__ cp)
{
  double T;

  {
    int const j = 0;

    T = fmin(10000.0, fmax(100.0, phi[0]));
    for (int k = 0; k <= 8; ++k)
    {
      if (T < T_mid[k])
        cp[k] = 8314.4621 * (T * (T * (T * (T * a_lo[4 + 7 * k] + a_lo[3 + 7 * k]) + a_lo[2 + 7 * k]) + a_lo[1 + 7 * k]) + a_lo[7 * k]);
      if (!(T < T_mid[k]))
        cp[k] = 8314.4621 * (T * (T * (T * (T * a_hi[4 + 7 * k] + a_hi[3 + 7 * k]) + a_hi[2 + 7 * k]) + a_hi[1 + 7 * k]) + a_hi[7 * k]);
    }
  }
}

void chem_utils(double const *__restrict__ t, double *__restrict__ b, double const *__restrict__ phi, double *__restrict__ h, double *__restrict__ cp, double *__restrict__ rwk)
{
    eval_b(phi, b);
    
    eval_h(phi, h);
    
    eval_cp(phi, cp);
}
