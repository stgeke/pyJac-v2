#include "jacobian_driver.hpp"
#include <math.h>
#include <stdlib.h>




void copy_in(double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ P_arr_local, double *__restrict__ phi_local, int const problem_size, int const driver_offset)
{
  {
    int const j = 0;

    if (j + driver_offset < problem_size)
      for (int i = 0; i <= 9; ++i)
      {
        phi_local[i] = phi[10 * (j + driver_offset) + i];
        P_arr_local[0] = P_arr[j + driver_offset];
      }
  }
}
void driver(int const problem_size, double const *__restrict__ t, double const *__restrict__ P_arr_local, double const *__restrict__ phi_local, double *__restrict__ jac_local, double *__restrict__ rwk)
{
  {
    int const i = 0;

    {
      int const j = 0;

      jacobian(t, P_arr_local, phi_local, jac_local, rwk);
    }
  }
}
void copy_out(double *__restrict__ jac, double const *__restrict__ jac_local, int const problem_size, int const driver_offset)
{
  {
    int const j = 0;

    if (j + driver_offset < problem_size)
      for (int i = 0; i <= 9; ++i)
        for (int i_0 = 0; i_0 <= 9; ++i_0)
          jac[100 * (j + driver_offset) + 10 * i + i_0] = jac_local[10 * i + i_0];
  }
}

void jacobian_driver(int const problem_size, double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ jac, double *__restrict__ rwk)
{
    
    #pragma omp parallel
    {
        // note: work_size unpacking _must_ be done in a parallel section in
        // order to get correct values from omp_get_num_threads()
        double* __restrict__ P_arr_local = rwk + 260 * work_size + 1 * omp_get_thread_num();
        double* __restrict__ phi_local = rwk + 261 * work_size + 10 * omp_get_thread_num();
        double* __restrict__ jac_local = rwk + 271 * work_size + 100 * omp_get_thread_num();
        double* __restrict__ t = NULL;
        #pragma omp for
        for (int driver_offset = 0; driver_offset < problem_size; driver_offset += 1)
        {
            copy_in(P_arr, phi, P_arr_local, phi_local, problem_size, driver_offset);
    
            driver(problem_size, t, P_arr_local, phi_local, jac_local, rwk);
    
            copy_out(jac, jac_local, problem_size, driver_offset);
        }
    }
}
