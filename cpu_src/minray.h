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

#define VERSION "0"

#define NONE 0
#define VACUUM 1
#define REFLECTIVE 2

#define SMALL 1
#define MEDIUM 2
#define LARGE 3

//#define BUMP 1.0e-11
//#define BUMP 1.0e-11
//#define BUMP 5.0e-6
#define BUMP 1.0e-6
//#define BUMP 1.0e-4

#define RT_FLOAT float

typedef struct{
  RT_FLOAT distance_to_surface;
  RT_FLOAT surface_normal_x;
  RT_FLOAT surface_normal_y;
  int last_surface;
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
  RT_FLOAT length_per_dimension;
  uint64_t n_rays;
  RT_FLOAT distance_per_ray;
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
  RT_FLOAT cell_expected_track_length;
  RT_FLOAT inverse_total_track_length;
  RT_FLOAT cell_width;
  RT_FLOAT inverse_length_per_dimension;
  RT_FLOAT inverse_cell_width;
  uint64_t n_cells;
  int n_iterations;
  int plotting_enabled;
  RT_FLOAT cell_volume;
  int validation_problem_id;
} Parameters;

typedef struct{
  float * angular_flux;
  RT_FLOAT * location_x;
  RT_FLOAT * location_y;
  RT_FLOAT * direction_x;
  RT_FLOAT * direction_y;
  int * cell_id;
  int * last_surface;
} RayData;

typedef struct{
  int * n_intersections;
  int * cell_ids;
  RT_FLOAT * distances;
  int * did_vacuum_reflects;
} IntersectionData;

typedef struct{
  float * isotropic_source;
  float * new_scalar_flux;
  float * old_scalar_flux;
  float * scalar_flux_accumulator;
  int   * hit_count;
  float * fission_rate;
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
  RT_FLOAT runtime_total;
  RT_FLOAT runtime_transport_sweep;
  RT_FLOAT k_eff;
  RT_FLOAT k_eff_std_dev;
} SimulationResult;

// io.c
Parameters read_CLI(int argc, char * argv[]);
ReadOnlyData load_2D_C5G7_XS(Parameters P);
void plot_3D_vtk(Parameters P, float * scalar_flux_accumulator, int * material_id);
void print_user_inputs(Parameters P);
int print_results(Parameters P, SimulationResult SR);
void print_status_data(int iter, RT_FLOAT k_eff, RT_FLOAT percent_missed, int is_active_region, RT_FLOAT k_eff_total_accumulator, RT_FLOAT k_eff_sum_of_squares_accumulator, int n_active_iterations);
void center_print(const char *s, int width);
void border_print(void);
void print_ray_tracing_buffer(Parameters P, SimulationData SD);
void print_ray(RT_FLOAT x, RT_FLOAT y, RT_FLOAT x_dir, RT_FLOAT y_dir, int cell_id);

// simulation.c
SimulationResult run_simulation(Parameters P, SimulationData SD);
void transport_sweep(Parameters P, SimulationData SD);
void update_isotropic_sources(Parameters P, SimulationData SD, RT_FLOAT k_eff);
void normalize_scalar_flux(Parameters P, SimulationData SD);
void add_source_to_scalar_flux(Parameters P, SimulationData SD);
void compute_cell_fission_rates(Parameters P, SimulationData SD, float * scalar_flux);
RT_FLOAT reduce_sum_float(float * a, int size);
int reduce_sum_int(int * a, int size);
RT_FLOAT compute_k_eff(Parameters P, SimulationData SD, RT_FLOAT old_k_eff);
RT_FLOAT check_hit_rate(int * hit_count, int n_cells);

// rand.c
RT_FLOAT LCG_random_double(uint64_t * seed);
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

// init.c
SimulationData initialize_simulation(Parameters P);
void initialize_rays(Parameters P, SimulationData SD);
void initialize_fluxes(Parameters P, SimulationData SD);
size_t estimate_memory_usage(Parameters P);

// utils.c
RT_FLOAT get_time(void);
void ptr_swap(float ** a, float ** b);
void compute_statistics(RT_FLOAT sum, RT_FLOAT sum_of_squares, int n, RT_FLOAT * sample_mean, RT_FLOAT * std_dev_of_sample_mean);
int validate_results(int validation_problem_id, RT_FLOAT k_eff);

// ray_trace_kernel.c
void ray_trace_kernel(Parameters P, SimulationData SD, RayData rayData, uint64_t ray_id);
TraceResult cartesian_ray_trace(RT_FLOAT x, RT_FLOAT y, RT_FLOAT cell_width, int x_idx, int y_idx, RT_FLOAT x_dir, RT_FLOAT y_dir, int last_surface);
CellLookup find_cell_id(Parameters P, RT_FLOAT x, RT_FLOAT y, int last_surface, int idx_x, int idx_y);

// Other kernel files
void update_isotropic_sources_kernel(Parameters P, SimulationData SD, int cell, int energy_group_in, RT_FLOAT inverse_k_eff);
void flux_attenuation_kernel(Parameters P, SimulationData SD, uint64_t ray_id, int energy_group);
void normalize_scalar_flux_kernel(Parameters P, float * new_scalar_flux, int cell, int energy_group);
void add_source_to_scalar_flux_kernel(Parameters P, SimulationData SD, int cell, int energy_group);
void compute_cell_fission_rates_kernel(Parameters P, SimulationData SD, float * scalar_flux, int cell);
