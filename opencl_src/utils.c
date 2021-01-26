#include "minray.h"
#include<sys/time.h>

double get_time(void)
{
  #ifdef OPENMP
  return omp_get_wtime();
  #endif

/*
  time_t time;
  time = clock();

  return (double) time / (double) CLOCKS_PER_SEC;
*/

  struct timeval start;
  gettimeofday(&start, NULL);

  double time_us = start.tv_sec * 1000000 + start.tv_usec;
  double time = time_us * 1e-6;

  return time;

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

int validate_results(int validation_problem_id, double k_eff)
{
  if(validation_problem_id)
  {
    double expected_results[3] = {0.31918, 1.19311, 1.18600};
    double k_eff_expected = expected_results[validation_problem_id - 1];
    double delta = fabs(k_eff - k_eff_expected);
    if( delta < 5.0e-5 )
    {
      printf("Validation Test                   = Passed\n");
      return 0;
    }
    else
    {
      printf("Validation Test                   = Failed\n");
      printf("    Expected k-eff                = %.7lf\n", k_eff_expected);
      printf("    Simulation k-eff              = %.7lf\n", k_eff);
      printf("    Delta                         = %.2le\n", delta);
      return 1;
    }
  }
  return 0;
}
