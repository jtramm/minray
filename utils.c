#include "minray.h"

double get_time(void)
{
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

void compute_statistics(double sum, double sum_of_squares, int n, double * sample_mean, double * std_dev_of_sample_mean)
{
  *std_dev_of_sample_mean = sqrt( (sum_of_squares - sum * sum / n ) / n);
  *std_dev_of_sample_mean /= sqrt(n);
  *sample_mean = sum / n;
}
