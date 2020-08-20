#include "minray.h"

void ray_trace_kernel(Parameters P, SimulationData SD)
{
  for( int r = 0; r < P.n_rays; r++ )
  {

  }
}

void transport_sweep(Parameters P, SimulationData SD)
{
  // Ray Trace Kernel
  ray_trace_kernel(P, SD);

  // Flux Attenuate Kernel

}

void run_simulation(Parameters P, SimulationData SD)
{
  double k-eff = 1.0;

  initialize_rays(P, SD);

  for( int iter = 0; iter < P.n_iterations; iter++ )
  {
    // Update Source
    
    // Set flux tally to zero

    // Transport Sweep

    // Normalize Scalar Flux

    // Add Source to Flux

    // Compute K-eff
  }
}

int main(int argc, char * argv[])
{
  Parameters P = read_CLI(argc, argv);
  SimulationData SD = initialize_simulation(P);

  return 0;
}
