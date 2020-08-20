#include "minray.h"

RayData initialize_ray_data(Parameters P)
{
  RayData rayData;

  size_t sz = P.n_rays * P.n_energy_groups * sizeof(float);
  rayData.angular_flux = (float *) malloc(sz);

  sz = P.n_rays * sizeof(double);
  rayData.location_x  = (double *) malloc(sz);
  rayData.location_y  = (double *) malloc(sz);
  rayData.location_z  = (double *) malloc(sz);
  rayData.direction_x = (double *) malloc(sz);
  rayData.direction_y = (double *) malloc(sz);
  rayData.direction_z = (double *) malloc(sz);
  
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
  CD.scalar_flux_accumulators = (float *) malloc(sz);

  sz = P.n_cells * sizeof(int);
  CD.hit_count = (int *) malloc(sz);

  return CD;
}

SimulationData initialize_simulation(Parameters P)
{
  ReadOnlyData ROD = load_2D_C5G7_XS(P);
  //TODO: Copy this to "device"
  
  ReadWriteData RWD;
  RWD.intersectionData = initialize_intersection_data(P);
  RWD.rayData          = initialize_ray_data(P);
  RWD.cellData         = initialize_cell_data(P);

  SimulationData SD;
  SD.readOnlyData  = ROD;
  SD.readWriteData = RWD;

  return SD;
}

#define PRNG_SAMPLES_PER_RAY 4
void initialize_ray_kernel(uint64_t base_seed, int ray_id, double length_per_dimension, RayData RD)
{
    uint64_t offset = ray_id * PRNG_SAMPLES_PER_RAY;
    uint64_t seed = fast_forward_LCG(base_seed, offset);
    RD.location_x[ray_id] = LCG_random_double(&seed) * length_per_dimension;
    RD.location_y[ray_id] = LCG_random_double(&seed) * length_per_dimension;
    RD.location_z[ray_id] = 0.0;

    // Sample Angle
    double theta = LCG_random_double(&seed) * 2.0 * M_PI;
    double z = -1.0 + 2.0 * LCG_random_double(&seed);
    double zo = sqrt(1.0 - z*z);

    // Spherical conversion
    RD.direction_x[ray_id] = zo * cosf(theta);
    RD.direction_y[ray_id] = zo * sinf(theta);
    RD.direction_z[ray_id] = z;

    // TODO: Normalize
}  

void initialize_rays(Parameters P, SimulationData SD)
{
  // Sample all rays in space and angle
  for( int r = 0; r < P.n_rays; r++ )
  {
    initialize_ray_kernel(P.seed, r, P.length_per_dimension, SD.readWriteData.rayData);
  }
}

void initialize_fluxes(Parameters P, SimulationData SD)
{
  // Old scalar flux set to 1.0
  for( int cell = 0; cell < P.n_cells; cell++ )
  {
    for( int g = 0; g < P.n_energy_groups; g++ ) 
    {
      SD.readWriteData.cellData.old_scalar_flux[         cell * P.n_energy_groups + g] = 1.0;
      SD.readWriteData.cellData.scalar_flux_accumulators[cell * P.n_energy_groups + g] = 1.0;
    }
  }
  
  // Set all starting angular fluxes to 0.0
  memset(SD.readWriteData.rayData.angular_flux, 0, P.n_rays * P.n_energy_groups * sizeof(float));
}
