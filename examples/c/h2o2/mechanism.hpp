#ifndef MECHANISM_HPP
#define MECHANISM_HPP
#ifdef _OPENMP
 #include <omp.h>
#else
 #warning 'OpenMP not found! Unexpected results may occur if using more than one thread.'
 #define omp_get_num_threads() (1)
#endif
#define NS (9)
#define NR (28)
#define NN (10)
#endif