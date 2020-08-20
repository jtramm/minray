#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include<string.h>

#define NONE 0
#define VACUUM 1
#define REFLECTIVE 2

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
  //double height;
  int n_rays;
  double distance_per_ray;
  int n_inactive_iterations;
  int n_active_iterations;
  uint64_t seed;
  int n_materials;
  int n_energy_groups;
  int max_intersections_per_ray;
  int boundary_conditions[5];
  int boundary_condition_x_positive;
  int boundary_condition_x_negative;
  int boundary_condition_y_positive;
  int boundary_condition_y_negative;
  // Derived
  double cell_expected_track_length;
  double cell_width;
  double inverse_length_per_dimension;
  double inverse_cell_width;
  int n_cells;
  int n_iterations;
} Parameters;

typedef struct{
  float * angular_flux;
  double * location_x;
  double * location_y;
  //double * location_z;
  double * direction_x;
  double * direction_y;
  //double * direction_z;
  int * cell_id;
} RayData;

typedef struct{
  int * n_intersections;
  int * cell_ids;
  double * distances;
  int * did_vacuum_reflects;
} IntersectionData;

typedef struct{
  float * isotropic_source;
  float * new_scalar_flux;
  float * old_scalar_flux;
  float * scalar_flux_accumulators;
  int   * hit_count;
} CellData;

typedef struct{
  CellData cellData;
  RayData rayData;
  IntersectionData intersectionData;
} ReadWriteData;

typedef struct{
  ReadOnlyData  readOnlyData;
  ReadWriteData readWriteData;
} SimulationData;

// io.c
Parameters read_CLI(int argc, char * argv[]);
ReadOnlyData load_2D_C5G7_XS(Parameters P);

// rand.c
double LCG_random_double(uint64_t * seed);
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

// init.c
SimulationData initialize_simulation(Parameters P);
void initialize_rays(Parameters P, SimulationData SD);
void initialize_fluxes(Parameters P, SimulationData SD);

// ray_trace_kernel.c
void ray_trace_kernel(Parameters P, SimulationData SD, RayData rayData, int ray_id);
double cartesian_ray_trace(double x, double y, double cell_width, int x_idx, int y_idx, double x_dir, double y_dir);
