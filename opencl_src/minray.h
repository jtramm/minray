#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<time.h>

//#define CL_TARGET_OPENCL_VERSION 200
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#define MAX_SOURCE_SIZE (0x100000)

#ifdef ALGORITHM_B
#include"neighbor_list_b.h"
#endif

#ifdef ALGORITHM_H
#include"neighbor_list_h.h"
#endif

#ifdef ALGORITHM_J
#include"neighbor_list_j.h"
#endif

#include "parameters.h"

#define VERSION "0"

#define SMALL 1
#define MEDIUM 2
#define LARGE 3

typedef struct{
  Node * nodes;
  int * idx;
  cl_mem d_nodes;
  size_t sz_nodes;
  cl_mem d_idx;
  size_t sz_idx;
} NodePool;
NodePool nl_init_nodePool(int n_cells);

typedef struct{
  cl_kernel ray_trace_kernel;
  cl_kernel flux_attenuation_kernel;          
  cl_kernel update_isotropic_sources_kernel;  
  cl_kernel normalize_scalar_flux_kernel;   
  cl_kernel add_source_to_scalar_flux_kernel; 
  cl_kernel compute_cell_fission_rates_kernel;
  cl_kernel reduce_int_kernel;
  cl_kernel reduce_float_kernel;
} Kernels;

typedef struct{
  cl_platform_id platform_id;
  cl_device_id device_id;
  cl_context context;
  cl_command_queue command_queue;
  Kernels kernels;
} OpenCLInfo;

typedef struct{
  int * material_id;
  float * nu_Sigma_f;
  float * Sigma_f;
  float * Sigma_t;
  float * Sigma_s;
  float * Chi;
  cl_mem d_material_id;
  cl_mem d_nu_Sigma_f;
  cl_mem d_Sigma_f;
  cl_mem d_Sigma_t;
  cl_mem d_Sigma_s;
  cl_mem d_Chi;
  size_t sz_material_id;
  size_t sz_nu_Sigma_f;
  size_t sz_Sigma_f;
  size_t sz_Sigma_t;
  size_t sz_Sigma_s;
  size_t sz_Chi;
} ReadOnlyData;


typedef struct{
  float  * angular_flux;
  double * location_x;
  double * location_y;
  double * direction_x;
  double * direction_y;
  int    * cell_id;
  cl_mem d_angular_flux;
  cl_mem d_location_x;
  cl_mem d_location_y;
  cl_mem d_direction_x;
  cl_mem d_direction_y;
  cl_mem d_cell_id;
  size_t sz_angular_flux;
  size_t sz_location_x;
  size_t sz_location_y;
  size_t sz_direction_x;
  size_t sz_direction_y;
  size_t sz_cell_id;
} RayData;

typedef struct{
  int    * n_intersections;
  int    * cell_ids;
  double * distances;
  int    * did_vacuum_reflects;
  cl_mem d_n_intersections;
  cl_mem d_cell_ids;
  cl_mem d_distances;
  cl_mem d_did_vacuum_reflects;
  size_t sz_n_intersections;
  size_t sz_cell_ids;
  size_t sz_distances;
  size_t sz_did_vacuum_reflects;
} IntersectionData;

typedef struct{
  float * isotropic_source;
  float * new_scalar_flux;
  float * old_scalar_flux;
  float * scalar_flux_accumulator;
  int   * hit_count;
  float * fission_rate;
  NeighborList * neighborList;
  cl_mem d_isotropic_source;
  cl_mem d_new_scalar_flux;
  cl_mem d_old_scalar_flux;
  cl_mem d_scalar_flux_accumulator;
  cl_mem d_hit_count;
  cl_mem d_fission_rate;
  cl_mem d_neighborList;
  size_t sz_isotropic_source;
  size_t sz_new_scalar_flux;
  size_t sz_old_scalar_flux;
  size_t sz_scalar_flux_accumulator;
  size_t sz_hit_count;
  size_t sz_fission_rate;
  size_t sz_neighborList;
} CellData;

typedef struct{
  CellData cellData;
  RayData rayData;
  IntersectionData intersectionData;
  NodePool nodePool;
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
SimulationResult run_simulation(OpenCLInfo * CL, Parameters P, SimulationData SD);
void transport_sweep(OpenCLInfo * CL, Parameters P, SimulationData SD);
void update_isotropic_sources(OpenCLInfo * CL, Parameters P, SimulationData SD, double k_eff);
void normalize_scalar_flux(OpenCLInfo * CL, Parameters P, SimulationData SD);
void add_source_to_scalar_flux(OpenCLInfo * CL, Parameters P, SimulationData SD);
double reduce_sum_float(float * a, int size);
int reduce_sum_int(int * a, int size);
double compute_k_eff(OpenCLInfo * CL, Parameters P, SimulationData SD, double old_k_eff);
void compute_cell_fission_rates(OpenCLInfo * CL, Parameters P, SimulationData SD, double utility_variable);
double check_hit_rate(OpenCLInfo * CL, SimulationData SD, int n_cells);
float reduce_fission_rates(OpenCLInfo *CL, SimulationData SD, int n_cells);
int reduce_intersections(OpenCLInfo *CL, SimulationData SD, int n_rays);

// rand.c
double LCG_random_double(uint64_t * seed);
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

// init.c
SimulationData initialize_simulation(Parameters P);
void initialize_rays(Parameters P, SimulationData SD);
void initialize_fluxes(Parameters P, SimulationData SD);
size_t estimate_memory_usage(Parameters P);
void initialize_device_data(SimulationData * SD, OpenCLInfo * CL);
void initialize_kernels(OpenCLInfo * CL);
void load_kernel_arguments(Parameters * P, SimulationData * SD, OpenCLInfo * CL);

// utils.c
double get_time(void);
void ptr_swap(float ** a, float ** b);
void compute_statistics(double sum, double sum_of_squares, int n, double * sample_mean, double * std_dev_of_sample_mean);
int validate_results(int validation_problem_id, double k_eff);

// cl_utils.c
const char *getErrorString(cl_int error);
void check(cl_int error);
void printCompilerError( cl_program program, cl_device_id device );
void print_single_info( cl_platform_id platform, cl_device_id device);
void print_opencl_info(void);
OpenCLInfo initialize_device(int user_platform_id, int user_device_id);
cl_mem copy_array_to_device(OpenCLInfo * CL, cl_mem_flags mem_flags, void * array, size_t sz);
cl_kernel compile_kernel(OpenCLInfo * CL, char * kernel_name);
void set_kernel_arguments(cl_kernel * kernel, int argc, size_t * arg_sz, void ** args);
void clear_array(OpenCLInfo * CL, cl_mem * buffer, size_t sz);
void copy_array_from_device(OpenCLInfo * CL, cl_mem * d_array, void * h_array, size_t sz);
