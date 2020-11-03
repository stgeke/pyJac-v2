#ifndef VECTORIZATION_HPP
#define VECTORIZATION_HPP
#ifdef _OPENMP
 #include <omp.h>
#else
 #warning 'OpenMP not found! Unexpected results may occur if using more than one thread.'
 #define omp_get_num_threads() (1)
#endif
#include "mechanism.hpp"
#endif