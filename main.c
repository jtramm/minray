#include "minray.h"

typedef struct{
  int * material_id;
  float * Sigma_f;
  float * Sigma_t;
  float * Sigma_s;
  float * Chi;
} ReadOnlyData;

typedef struct{
  // User Inputs
  int n_cells_per_dimension;
  double length_per_dimension;
  double height;
  int n_rays;
  double distance_per_ray;
  int n_inactive_iterations;
  int n_active_iterations;
  uint64_t seed;
  int n_materials;
  int n_energy_groups;
  // Derived
  double cell_volume;
  double cell_expected_track_length;
  double cell_width;
  int n_cells;
} Parameters;

typedef struct{
  float * delta_psi_tally;
  float * isotropic_source;
  float * scalar_flux;
  float * scalar_flux_accumulators;
  int   * hit_count;
} ReadWriteData;

Parameters read_CLI(int argc, char * argv[])
{
  Parameters P;
  // User Inputs
  P.n_cells_per_dimension = 100;
  P.length_per_dimension = 20.0;
  P.height = 100.0;
  P.n_rays = 5;
  P.distance_per_ray = 10.0;
  P.n_inactive_iterations = 10;
  P.n_active_iterations = 10;
  P.seed = 1337;
  P.n_materials = 5;
  P.n_energy_groups = 7;

  // Derived Values
  P.cell_width = P.length_per_dimension / P.n_cells_per_dimension;
  P.cell_volume = P.cell_width * P.cell_width * P.height;
  P.n_cells = n_cells_per_dimension * n_cells_per_dimension;
  P.cell_expected_track_length = (distance_per_ray * n_rays) / n_cells;
  
  return P;
}



int main(int argc, char * argv[])
{
  Parameters P = read_CLI(argc, argv);

  return 0;
}
