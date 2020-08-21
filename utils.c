#include "minray.h"

double get_time(void)
{
  #ifdef _MPI
  return MPI_Wtime();
  #endif

  #ifdef OPENMP
  return omp_get_wtime();
  #endif

  time_t time;
  time = clock();

  return (double) time / (double) CLOCKS_PER_SEC;
}

void ptr_swap(float ** a, float ** b)
{
  float * tmp = *a;
  *a = *b;
  *b = tmp;
}
