#ifndef JACOBIAN_DRIVER_HPP
#define JACOBIAN_DRIVER_HPP
#ifdef _OPENMP
 #include <omp.h>
#else
 #warning 'OpenMP not found! Unexpected results may occur if using more than one thread.'
 #define omp_get_num_threads() (1)
#endif
#include "mechanism.hpp"
#include "jacobian.hpp"

#ifndef work_size
    #define work_size (omp_get_num_threads())
#endif

void copy_in(double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ P_arr_local, double *__restrict__ phi_local, int const problem_size, int const driver_offset);
void driver(int const problem_size, double const *__restrict__ t, double const *__restrict__ P_arr_local, double const *__restrict__ phi_local, double *__restrict__ jac_local, double *__restrict__ rwk);
void copy_out(double *__restrict__ jac, double const *__restrict__ jac_local, int const problem_size, int const driver_offset);
void jacobian_driver(int const problem_size, double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ jac, double *__restrict__ rwk);
#endif