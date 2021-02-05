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
  sz += (P.n_cells * P.n_energy_groups * sizeof(float))*4;
  sz += P.n_cells * sizeof(float);
  sz += P.n_cells * sizeof(int);
  sz += P.n_cells * sizeof(NeighborList);
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

  sz = P.n_rays * sizeof(double);
  rayData.location_x  = (double *) malloc(sz);
  rayData.location_y  = (double *) malloc(sz);
  rayData.direction_x = (double *) malloc(sz);
  rayData.direction_y = (double *) malloc(sz);
  
  sz = P.n_rays * sizeof(int);
  rayData.cell_id  = (int *) malloc(sz);

  return rayData;
}

IntersectionData initialize_intersection_data(Parameters P)
{
  IntersectionData intersectionData;

  size_t sz = P.n_rays * P.max_intersections_per_ray * sizeof(int);
  intersectionData.n_intersections     = (int *) malloc(sz);
  intersectionData.cell_ids            = (int *) malloc(sz);
  intersectionData.did_vacuum_reflects = (int *) malloc(sz);

  sz        = P.n_rays * P.max_intersections_per_ray * sizeof(double);
  intersectionData.distances = (double *) malloc(sz);

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

  sz = P.n_cells * sizeof(float);
  CD.fission_rate             = (float *) malloc(sz);

  sz = P.n_cells * sizeof(int);
  CD.hit_count = (int *) malloc(sz);

  sz = P.n_cells * sizeof(NeighborList);
  CD.neighborList = (NeighborList *) malloc(sz);

  for( int i = 0; i < P.n_cells; i++ )
    nl_init(CD.neighborList + i);

  return CD;
}

SimulationData initialize_simulation(Parameters P)
{
  border_print();
  center_print("INITIALIZATION", 79);
  border_print();

  printf("Initializing read only data...\n");
  ReadOnlyData ROD = load_2D_C5G7_XS(P);
  
  printf("Initializing read/write data...\n");
  ReadWriteData RWD;
  RWD.intersectionData = initialize_intersection_data(P);
  RWD.rayData          = initialize_ray_data(P);
  RWD.cellData         = initialize_cell_data(P);
  RWD.neighborListPool = nl_init_pool(P.n_cells);

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
