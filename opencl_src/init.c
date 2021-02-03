#include "minray.h"

size_t estimate_memory_usage(Parameters P)
{
  size_t sz = 0;
  // Ray Data
  sz += P.n_rays * P.n_energy_groups * sizeof(float);
  sz += (P.n_rays * sizeof(double)) * 4;
  sz += P.n_rays * sizeof(int);
  // Intersection Data
  sz += (P.n_rays * P.max_intersections_per_ray * sizeof(int))*3;
  sz += P.n_rays * P.max_intersections_per_ray * sizeof(double);
  // Cell Data
  sz += (P.n_cells * P.n_energy_groups * sizeof(float))*5;
  sz += P.n_cells * sizeof(int);
  // XS Data
  sz += P.n_materials * P.n_energy_groups * sizeof(float)*4;
  sz += P.n_materials * P.n_energy_groups * P.n_energy_groups * sizeof(float);
  sz += P.n_cells * sizeof(int);
  return sz;
}

RayData initialize_ray_data(Parameters P)
{
  RayData rayData;

  size_t sz = P.n_rays * P.n_energy_groups * sizeof(float);
  rayData.angular_flux = (float *) malloc(sz);
  rayData.sz_angular_flux = sz;

  sz = P.n_rays * sizeof(double);
  rayData.location_x  = (double *) malloc(sz);
  rayData.location_y  = (double *) malloc(sz);
  rayData.direction_x = (double *) malloc(sz);
  rayData.direction_y = (double *) malloc(sz);
  rayData.sz_location_x  = sz;
  rayData.sz_location_y  = sz;
  rayData.sz_direction_x = sz;
  rayData.sz_direction_y = sz;
  
  sz = P.n_rays * sizeof(int);
  rayData.cell_id  = (int *) malloc(sz);
  rayData.sz_cell_id  = sz;

  return rayData;
}

IntersectionData initialize_intersection_data(Parameters P)
{
  IntersectionData intersectionData;

  size_t sz = P.n_rays * P.max_intersections_per_ray * sizeof(int);
  intersectionData.n_intersections     = (int *) malloc(sz);
  intersectionData.cell_ids            = (int *) malloc(sz);
  intersectionData.did_vacuum_reflects = (int *) malloc(sz);
  intersectionData.sz_n_intersections     = sz;
  intersectionData.sz_cell_ids            = sz;
  intersectionData.sz_did_vacuum_reflects = sz;

  sz        = P.n_rays * P.max_intersections_per_ray * sizeof(double);
  intersectionData.distances = (double *) malloc(sz);
  intersectionData.sz_distances = sz;

  return intersectionData;
}

CellData initialize_cell_data(Parameters P)
{
  CellData CD;

  size_t sz = P.n_cells * P.n_energy_groups * sizeof(float);
  CD.isotropic_source         = (float *) malloc(sz);
  CD.new_scalar_flux          = (float *) malloc(sz);
  CD.old_scalar_flux          = (float *) malloc(sz);
  CD.scalar_flux_accumulator  = (float *) malloc(sz);
  CD.sz_isotropic_source         = sz;
  CD.sz_new_scalar_flux          = sz;
  CD.sz_old_scalar_flux          = sz;
  CD.sz_scalar_flux_accumulator  = sz;

  sz = P.n_cells * sizeof(float);
  CD.fission_rate             = (float *) malloc(sz);
  CD.sz_fission_rate             = sz;

  sz = P.n_cells * sizeof(int);
  CD.hit_count = (int *) malloc(sz);
  CD.sz_hit_count = sz;
  
  sz = P.n_cells * sizeof(NeighborList);
  CD.neighborList = (NeighborList *) malloc(sz);

  for( int i = 0; i < P.n_cells; i++ )
    nl_init(CD.neighborList + i);

  CD.sz_neighborList = sz;

  return CD;
}

SimulationData initialize_simulation(Parameters P)
{
  border_print();
  center_print("HOST INITIALIZATION", 79);
  border_print();

  printf("Initializing read only data...\n");
  ReadOnlyData ROD = load_2D_C5G7_XS(P);
  
  printf("Initializing read/write data...\n");
  ReadWriteData RWD;
  RWD.intersectionData = initialize_intersection_data(P);
  RWD.rayData          = initialize_ray_data(P);
  RWD.cellData         = initialize_cell_data(P);
  RWD.sz_vectorPool    = P.n_cells * sizeof(int) * AVG_NEIGHBORS_PER_CELL; 
  RWD.vectorPool       = (int *) malloc(RWD.sz_vectorPool);
  for( int i = 0; i < P.n_cells * AVG_NEIGHBORS_PER_CELL; i++ )
    RWD.vectorPool[i] = -1;
  RWD.sz_vectorPool_idx = sizeof(int);
  RWD.vectorPool_idx    = (int *) malloc(RWD.sz_vectorPool_idx);
  *RWD.vectorPool_idx = 0;

  SimulationData SD;
  SD.readOnlyData  = ROD;
  SD.readWriteData = RWD;
  
  border_print();

  return SD;
}

#define PRNG_SAMPLES_PER_RAY 10
void initialize_ray_kernel(uint64_t base_seed, int ray_id, double length_per_dimension, int n_cells_per_dimension, double inverse_cell_width, RayData RD)
{
    uint64_t offset = ray_id * PRNG_SAMPLES_PER_RAY;
    uint64_t seed = fast_forward_LCG(base_seed, offset);
    RD.location_x[ray_id] = LCG_random_double(&seed) * length_per_dimension;
    RD.location_y[ray_id] = LCG_random_double(&seed) * length_per_dimension;

    // Sample azimuthal angle
    double theta = LCG_random_double(&seed) * 2.0 * M_PI;

    // Sample polar angle
    double z = -1.0 + 2.0 * LCG_random_double(&seed);

    // If polar angle approaches unity (i.e., very steep), this can cause numerical instability.
    // To fix this, for polar angles ~1.0 we will just resample.
    while(fabs(z) > 0.9999)
      z = -1.0 + 2.0 * LCG_random_double(&seed);

    // Spherical conversion
    double zo = sqrt(1.0 - z*z);
    double x = zo * cosf(theta);
    double y = zo * sinf(theta);

    // Normalize Direction
    double inverse = 1.0 / sqrt( x*x + y*y + z*z );
    x *= inverse;
    y *= inverse;
    z *= inverse;
    
    // Compute Starting Cell ID
    int x_idx = RD.location_x[ray_id] * inverse_cell_width;
    int y_idx = RD.location_y[ray_id] * inverse_cell_width;
    int cell_id = y_idx * n_cells_per_dimension + x_idx;

    // Store sampled ray data
    RD.cell_id[    ray_id] = cell_id; 
    RD.direction_x[ray_id] = x;
    RD.direction_y[ray_id] = y;
}  

void initialize_rays(Parameters P, SimulationData SD)
{
  // Sample all rays in space and angle
  for( int r = 0; r < P.n_rays; r++ )
  {
    initialize_ray_kernel(P.seed, r, P.length_per_dimension, P.n_cells_per_dimension, P.inverse_cell_width, SD.readWriteData.rayData);
  }
}

void initialize_fluxes(Parameters P, SimulationData SD)
{
  // Old scalar fluxes set to 1.0
  for( int cell = 0; cell < P.n_cells; cell++ )
  {
    for( int g = 0; g < P.n_energy_groups; g++ ) 
    {
      SD.readWriteData.cellData.old_scalar_flux[cell * P.n_energy_groups + g] = 1.0;
    }
  }

  // Set scalar flux accumulators to 0.0
  memset(SD.readWriteData.cellData.scalar_flux_accumulator, 0, P.n_cells * P.n_energy_groups * sizeof(float));

  // Set all starting angular fluxes to 0.0
  memset(SD.readWriteData.rayData.angular_flux, 0, P.n_rays * P.n_energy_groups * sizeof(float));
}

void initialize_device_data(SimulationData * SD, OpenCLInfo * CL)
{
  // Copy Read Only Data
  cl_mem_flags mem_type = CL_MEM_READ_ONLY;
  SD->readOnlyData.d_material_id = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.material_id, SD->readOnlyData.sz_material_id);
  SD->readOnlyData.d_nu_Sigma_f  = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.nu_Sigma_f,  SD->readOnlyData.sz_nu_Sigma_f);
  SD->readOnlyData.d_Sigma_f     = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Sigma_f,     SD->readOnlyData.sz_Sigma_f);
  SD->readOnlyData.d_Sigma_t     = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Sigma_t,     SD->readOnlyData.sz_Sigma_t);
  SD->readOnlyData.d_Sigma_s     = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Sigma_s,     SD->readOnlyData.sz_Sigma_s);
  SD->readOnlyData.d_Chi         = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Chi,         SD->readOnlyData.sz_Chi);

  // Copy Read Write Data
  mem_type = CL_MEM_READ_WRITE;
  SD->readWriteData.rayData.d_angular_flux = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.angular_flux, SD->readWriteData.rayData.sz_angular_flux);
  SD->readWriteData.rayData.d_location_x   = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.location_x,   SD->readWriteData.rayData.sz_location_x);
  SD->readWriteData.rayData.d_location_y   = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.location_y,   SD->readWriteData.rayData.sz_location_y);
  SD->readWriteData.rayData.d_direction_x  = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.direction_x,  SD->readWriteData.rayData.sz_direction_x);
  SD->readWriteData.rayData.d_direction_y  = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.direction_y,  SD->readWriteData.rayData.sz_direction_y);
  SD->readWriteData.rayData.d_cell_id      = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.cell_id,      SD->readWriteData.rayData.sz_cell_id);
  
  SD->readWriteData.intersectionData.d_n_intersections     = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.n_intersections,     SD->readWriteData.intersectionData.sz_n_intersections);
  SD->readWriteData.intersectionData.d_cell_ids            = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.cell_ids,            SD->readWriteData.intersectionData.sz_cell_ids);
  SD->readWriteData.intersectionData.d_distances           = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.distances,           SD->readWriteData.intersectionData.sz_distances);
  SD->readWriteData.intersectionData.d_did_vacuum_reflects = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.did_vacuum_reflects, SD->readWriteData.intersectionData.sz_did_vacuum_reflects);
  
  SD->readWriteData.cellData.d_isotropic_source        = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.isotropic_source,        SD->readWriteData.cellData.sz_isotropic_source);
  SD->readWriteData.cellData.d_new_scalar_flux         = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.new_scalar_flux,         SD->readWriteData.cellData.sz_new_scalar_flux);
  SD->readWriteData.cellData.d_old_scalar_flux         = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.old_scalar_flux,         SD->readWriteData.cellData.sz_old_scalar_flux);
  SD->readWriteData.cellData.d_scalar_flux_accumulator = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.scalar_flux_accumulator, SD->readWriteData.cellData.sz_scalar_flux_accumulator);
  SD->readWriteData.cellData.d_hit_count               = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.hit_count,               SD->readWriteData.cellData.sz_hit_count);
  SD->readWriteData.cellData.d_fission_rate            = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.fission_rate,            SD->readWriteData.cellData.sz_fission_rate);
  SD->readWriteData.cellData.d_neighborList            = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.neighborList,            SD->readWriteData.cellData.sz_neighborList);

  // Copy vector pool data
  SD->readWriteData.d_vectorPool          = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.vectorPool,            SD->readWriteData.sz_vectorPool);
}

void initialize_kernels(OpenCLInfo * CL)
{
  CL->kernels.normalize_scalar_flux_kernel      = compile_kernel(CL, "normalize_scalar_flux_kernel");
  CL->kernels.ray_trace_kernel                  = compile_kernel(CL, "ray_trace_kernel");
  CL->kernels.flux_attenuation_kernel           = compile_kernel(CL, "flux_attenuation_kernel");
  CL->kernels.update_isotropic_sources_kernel   = compile_kernel(CL, "update_isotropic_sources_kernel");
  CL->kernels.add_source_to_scalar_flux_kernel  = compile_kernel(CL, "add_source_to_scalar_flux_kernel");
  CL->kernels.compute_cell_fission_rates_kernel = compile_kernel(CL, "compute_cell_fission_rates_kernel");
  CL->kernels.reduce_int_kernel                 = compile_kernel(CL, "reduce_int_kernel");
  CL->kernels.reduce_float_kernel               = compile_kernel(CL, "reduce_float_kernel");
}

#define NUM_ARGS 27
void load_kernel_arguments(Parameters * P, SimulationData * SD, OpenCLInfo * CL)
{
  printf("Loading static kernel arguments...\n");
  int argc;
  size_t arg_sz[NUM_ARGS];
  void * args[NUM_ARGS];

  argc = NUM_ARGS;

  // Set argument sizes
  arg_sz[0] = sizeof(double);
  arg_sz[1] = sizeof(Parameters);
  for( int i = 2; i < argc; i++ )
    arg_sz[i] = sizeof(cl_mem);

  double utility_variable = 0.0;

  // Set arguments
  args[ 0] = (void *) &utility_variable; 
  args[ 1] = (void *) P;

  args[ 2] = (void *) &SD->readOnlyData.d_material_id;
  args[ 3] = (void *) &SD->readOnlyData.d_nu_Sigma_f;
  args[ 4] = (void *) &SD->readOnlyData.d_Sigma_f;
  args[ 5] = (void *) &SD->readOnlyData.d_Sigma_t;
  args[ 6] = (void *) &SD->readOnlyData.d_Sigma_s;
  args[ 7] = (void *) &SD->readOnlyData.d_Chi;

  args[ 8] = (void *) &SD->readWriteData.rayData.d_angular_flux;
  args[ 9] = (void *) &SD->readWriteData.rayData.d_location_x;
  args[10] = (void *) &SD->readWriteData.rayData.d_location_y;
  args[11] = (void *) &SD->readWriteData.rayData.d_direction_x;
  args[12] = (void *) &SD->readWriteData.rayData.d_direction_y;
  args[13] = (void *) &SD->readWriteData.rayData.d_cell_id;

  args[14] = (void *) &SD->readWriteData.intersectionData.d_n_intersections;
  args[15] = (void *) &SD->readWriteData.intersectionData.d_cell_ids;
  args[16] = (void *) &SD->readWriteData.intersectionData.d_distances;
  args[17] = (void *) &SD->readWriteData.intersectionData.d_did_vacuum_reflects;

  args[18] = (void *) &SD->readWriteData.cellData.d_isotropic_source;
  args[19] = (void *) &SD->readWriteData.cellData.d_new_scalar_flux;
  args[20] = (void *) &SD->readWriteData.cellData.d_old_scalar_flux;
  args[21] = (void *) &SD->readWriteData.cellData.d_scalar_flux_accumulator;
  args[22] = (void *) &SD->readWriteData.cellData.d_hit_count;
  args[23] = (void *) &SD->readWriteData.cellData.d_fission_rate;
  args[24] = (void *) &SD->readWriteData.cellData.d_neighborList;
  args[25] = (void *) &SD->readWriteData.d_vectorPool;
  args[26] = (void *) &SD->readWriteData.d_vectorPool_idx;

  set_kernel_arguments(&CL->kernels.normalize_scalar_flux_kernel     , argc, arg_sz, args);
  set_kernel_arguments(&CL->kernels.ray_trace_kernel                 , argc, arg_sz, args); 
  set_kernel_arguments(&CL->kernels.flux_attenuation_kernel          , argc, arg_sz, args); 
  set_kernel_arguments(&CL->kernels.update_isotropic_sources_kernel  , argc, arg_sz, args); 
  set_kernel_arguments(&CL->kernels.add_source_to_scalar_flux_kernel , argc, arg_sz, args); 
  set_kernel_arguments(&CL->kernels.compute_cell_fission_rates_kernel, argc, arg_sz, args); 
}
