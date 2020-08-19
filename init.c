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
  intersectionData.n_intersections    = (int *) malloc(sz);
  intersectionData.cell_ids            = (int *) malloc(sz);
  intersectionData.did_vacuum_reflects = (int *) malloc(sz);

  sz        = P.n_rays * P.max_intersections_per_ray * sizeof(double);
  intersectionData.distances = (double *) malloc(sz);

  return intersectionData;
}

CellData initialize_cell_data(Parameters P)
{
  CellData CD;

  size_t sz = P.n_cells * sizeof(float);
  CD.delta_psi_tally          = (float *) malloc(sz);
  CD.isotropic_source         = (float *) malloc(sz);
  CD.scalar_flux              = (float *) malloc(sz);
  CD.scalar_flux_accumulators = (float *) malloc(sz);

  sz = P.n_cells * sizeof(int);
  CD.hit_count = (int *) malloc(sz);

  return CD;
}

SimulationData initialize_simulation(Parameters P)
{
  ReadOnlyData ROD = load_2D_C5G7_XS(P);
  
  ReadWriteData RWD;
  RWD.intersectionData = initialize_intersection_data(P);
  RWD.rayData          = initialize_ray_data(P);
  RWD.cellData         = initialize_cell_data(P);

  SimulationData SD;
  SD.readOnlyData  = ROD;
  SD.readWriteData = RWD;

  return SD;
}
