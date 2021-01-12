#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<time.h>
#ifdef OPENMP
#include<omp.h>
#endif

#ifdef ALGORITHM_A
#include"neighbor_list_a.h"
#endif

#ifdef ALGORITHM_B
#include"neighbor_list_b.h"
#endif

#ifdef ALGORITHM_C
#include"neighbor_list_c.h"
#endif

#ifdef ALGORITHM_D
#include"neighbor_list_d.h"
#endif

#ifdef ALGORITHM_E
#include"neighbor_list_e.h"
#endif

#ifdef ALGORITHM_F
#include"neighbor_list_f.h"
#endif

#define VERSION "0"

#define NONE 0
#define VACUUM 1
#define REFLECTIVE 2

#define SMALL 1
#define MEDIUM 2
#define LARGE 3

#define BUMP 1.0e-11

typedef struct{
  double distance_to_surface;
  double surface_normal_x;
  double surface_normal_y;
} TraceResult;

typedef struct{
  int cell_id;
  int cartesian_cell_idx_x;
  int cartesian_cell_idx_y;
  int boundary_condition;
} CellLookup;

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
  uint64_t n_rays;
  double distance_per_ray;
  int n_inactive_iterations;
  int n_active_iterations;
  uint64_t seed;
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
  uint64_t n_cells;
  int n_iterations;
  int plotting_enabled;
  double cell_volume;
  int validation_problem_id;
} Parameters;

typedef struct{
  float * angular_flux;
  double * location_x;
  double * location_y;
  double * direction_x;
  double * direction_y;
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
  float * scalar_flux_accumulator;
  int   * hit_count;
  float * fission_rate;
  NeighborList * neighborList;
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


typedef struct{
  uint64_t n_geometric_intersections;
  double runtime_total;
  double runtime_transport_sweep;
  double k_eff;
  double k_eff_std_dev;
} SimulationResult;

// io.c
Parameters read_CLI(int argc, char * argv[]);
ReadOnlyData load_2D_C5G7_XS(Parameters P);
void plot_3D_vtk(Parameters P, float * scalar_flux_accumulator, int * material_id);
void print_user_inputs(Parameters P);
int print_results(Parameters P, SimulationResult SR);
void print_status_data(int iter, double k_eff, double percent_missed, int is_active_region, double k_eff_total_accumulator, double k_eff_sum_of_squares_accumulator, int n_active_iterations);
void center_print(const char *s, int width);
void border_print(void);
void print_ray_tracing_buffer(Parameters P, SimulationData SD);
void print_ray(double x, double y, double x_dir, double y_dir, int cell_id);

// simulation.c
SimulationResult run_simulation(Parameters P, SimulationData SD);
void transport_sweep(Parameters P, SimulationData SD);
void update_isotropic_sources(Parameters P, SimulationData SD, double k_eff);
void normalize_scalar_flux(Parameters P, SimulationData SD);
void add_source_to_scalar_flux(Parameters P, SimulationData SD);
void compute_cell_fission_rates(Parameters P, SimulationData SD, float * scalar_flux);
double reduce_sum_float(float * a, int size);
int reduce_sum_int(int * a, int size);
double compute_k_eff(Parameters P, SimulationData SD, double old_k_eff);
double check_hit_rate(int * hit_count, int n_cells);

// rand.c
double LCG_random_double(uint64_t * seed);
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

// init.c
SimulationData initialize_simulation(Parameters P);
void initialize_rays(Parameters P, SimulationData SD);
void initialize_fluxes(Parameters P, SimulationData SD);
size_t estimate_memory_usage(Parameters P);

// utils.c
double get_time(void);
void ptr_swap(float ** a, float ** b);
void compute_statistics(double sum, double sum_of_squares, int n, double * sample_mean, double * std_dev_of_sample_mean);
int validate_results(int validation_problem_id, double k_eff);

// ray_trace_kernel.c
void ray_trace_kernel(Parameters P, SimulationData SD, RayData rayData, uint64_t ray_id);
CellLookup find_cell_id(Parameters P, double x, double y);
TraceResult cartesian_ray_trace(double x, double y, double cell_width, int x_idx, int y_idx, double x_dir, double y_dir);
CellLookup find_cell_id_using_neighbor_list(Parameters P, NeighborList * neighborList, double x, double y);

// Other kernel files
void update_isotropic_sources_kernel(Parameters P, SimulationData SD, int cell, int energy_group_in, double inverse_k_eff);
void flux_attenuation_kernel(Parameters P, SimulationData SD, uint64_t ray_id, int energy_group);
void normalize_scalar_flux_kernel(Parameters P, float * new_scalar_flux, int cell, int energy_group);
void add_source_to_scalar_flux_kernel(Parameters P, SimulationData SD, int cell, int energy_group);
void compute_cell_fission_rates_kernel(Parameters P, SimulationData SD, float * scalar_flux, int cell);
