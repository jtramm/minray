#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>

typedef struct{
  int * material_id;
  float * nu_Sigma_f;
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

Parameters read_CLI(int argc, char * argv[]);
ReadOnlyData load_2D_C5G7_XS(Parameters P);
