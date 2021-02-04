#define BUMP 1.0e-11
#define NONE 0
#define VACUUM 1
#define REFLECTIVE 2

#ifdef ALGORITHM_B
#include"neighbor_list_b.h"
#define ALGORITHM "-DALGORITHM_B"
#endif

#ifdef ALGORITHM_H
#include"neighbor_list_h.h"
#define ALGORITHM "-DALGORITHM_H"
#endif

#ifdef ALGORITHM_J
#include"neighbor_list_j.h"
#define ALGORITHM "-DALGORITHM_J -I."
#endif

#define ARGUMENTS double utility_variable, Parameters P, __global int * material_id, __global  float * nu_Sigma_f, __global  float * Sigma_f, __global  float * Sigma_t, __global  float * Sigma_s, __global  float * Chi, __global  float  * angular_flux, __global  double * location_x, __global  double * location_y, __global  double * direction_x, __global  double * direction_y, __global  int    * cell_id, __global  int    * n_intersections, __global  int    * cell_ids, __global  double * distances, __global  int    * did_vacuum_reflects, __global  float * isotropic_source, __global  float * new_scalar_flux, __global  float * old_scalar_flux, __global  float * scalar_flux_accumulator, __global  int   * hit_count, __global  float * fission_rate, __global NeighborList * neighborList, __global Node * nodePool_nodes, __global int * nodePool_idx

typedef struct{
  // User Inputs
  int n_cells_per_dimension;
  double length_per_dimension;
  ulong n_rays;
  double distance_per_ray;
  int n_inactive_iterations;
  int n_active_iterations;
  ulong seed;
  int n_materials;
  int n_energy_groups;
  int max_intersections_per_ray;
  int boundary_conditions[3][3];
  int boundary_condition_x_positive;
  int boundary_condition_x_negative;
  int boundary_condition_y_positive;
  int boundary_condition_y_negative;
  // Derived
  double cell_expected_track_length;
  double inverse_total_track_length;
  double cell_width;
  double inverse_length_per_dimension;
  double inverse_cell_width;
  ulong n_cells;
  int n_iterations;
  int plotting_enabled;
  double cell_volume;
  int validation_problem_id;
  int platform_id;
  int device_id;
  int n_nodes;
} Parameters;
