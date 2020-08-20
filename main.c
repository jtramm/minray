#include "minray.h"

void transport_sweep(Parameters P, SimulationData SD)
{
  // Ray Trace Kernel
  for( int r = 0; r < P.n_rays; r++ )
    ray_trace_kernel(P, SD, SD.readWriteData.rayData, r);

  // Flux Attenuate Kernel
}

void run_simulation(Parameters P, SimulationData SD)
{
  // k is the multiplication factor (or eigenvalue) we are trying to solve for.
  // The eigenvector is the scalar flux vector
  double k = 1.0;

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
  initialize_rays(P, SD);
  initialize_fluxes(P, SD);

  {
  double cell_width = P.cell_width;
  printf("cell_width = %.7lf\n", cell_width);
  int x_idx = 4;
  int y_idx = 6;
  double x = x_idx * cell_width + 0.001;
  double y = y_idx * cell_width + 0.001;
  double x_dir = sqrt(2.0) / 2.0;
  double y_dir = sqrt(2.0) / 2.0;
  double z_dir = 0;

  double dist = cartesian_ray_trace(x, y, cell_width, x_idx, y_idx, x_dir, y_dir, z_dir);
  printf("distance = %.5lf\n", dist);
  }
  {
  double cell_width = P.cell_width;
  printf("cell_width = %.7lf\n", cell_width);
  int x_idx = 4;
  int y_idx = 6;
  double x = x_idx * cell_width + 0.9999 * cell_width;
  double y = y_idx * cell_width + 0.9999 * cell_width;
  double x_dir = -sqrt(2.0) / 2.0;
  double y_dir = -sqrt(2.0) / 2.0;
  double z_dir = 0;

  double dist = cartesian_ray_trace(x, y, cell_width, x_idx, y_idx, x_dir, y_dir, z_dir);
  printf("distance = %.5lf\n", dist);
  }
  
  {
  double cell_width = P.cell_width;
  printf("cell_width = %.7lf\n", cell_width);
  int x_idx = 4;
  int y_idx = 6;
  double x = x_idx * cell_width + cell_width * 0.5;
  double y = y_idx * cell_width + 0.001;
  double x_dir = sqrt(2.0) / 2.0;
  double y_dir = sqrt(2.0) / 2.0;
  double z_dir = 0;

  double dist = cartesian_ray_trace(x, y, cell_width, x_idx, y_idx, x_dir, y_dir, z_dir);
  printf("distance = %.5lf\n", dist);
  }
  
  {
  double cell_width = P.cell_width;
  printf("cell_width = %.7lf\n", cell_width);
  int x_idx = 4;
  int y_idx = 6;
  double x = x_idx * cell_width + cell_width * 0.5;
  double y = y_idx * cell_width + 0.001;
  double x_dir = 0.00001;
  double y_dir = 0.99999;
  double z_dir = 0;

  double dist = cartesian_ray_trace(x, y, cell_width, x_idx, y_idx, x_dir, y_dir, z_dir);
  printf("distance = %.5lf\n", dist);
  }


  return 0;
}
