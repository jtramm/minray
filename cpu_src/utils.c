#include "minray.h"

RT_FLOAT get_time(void)
{
  #ifdef OPENMP
  return omp_get_wtime();
  #endif

  time_t time;
  time = clock();

  return (RT_FLOAT) time / (RT_FLOAT) CLOCKS_PER_SEC;
}

void ptr_swap(float ** a, float ** b)
{
  float * tmp = *a;
  *a = *b;
  *b = tmp;
}

void compute_statistics(RT_FLOAT sum, RT_FLOAT sum_of_squares, int n, RT_FLOAT * sample_mean, RT_FLOAT * std_dev_of_sample_mean)
{
  *std_dev_of_sample_mean = sqrt( (sum_of_squares - sum * sum / n ) / n);
  *std_dev_of_sample_mean /= sqrt(n);
  *sample_mean = sum / n;
}

int validate_results(int validation_problem_id, RT_FLOAT k_eff)
{
  if(validation_problem_id)
  {
    RT_FLOAT expected_results[3] = {0.31918, 1.19311, 1.18600};
    RT_FLOAT k_eff_expected = expected_results[validation_problem_id - 1];
    RT_FLOAT delta = fabs(k_eff - k_eff_expected);
    if( delta < 1.0e-5 )
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
