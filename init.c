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

  sz = P.n_rays * sizeof(double);
  rayData.location_x  = (double *) malloc(sz);
  rayData.location_y  = (double *) malloc(sz);
  //rayData.location_z  = (double *) malloc(sz);
  rayData.direction_x = (double *) malloc(sz);
  rayData.direction_y = (double *) malloc(sz);
 // rayData.direction_z = (double *) malloc(sz);
  
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
  CD.fission_rate             = (float *) malloc(sz);

  sz = P.n_cells * sizeof(int);
  CD.hit_count = (int *) malloc(sz);

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

  SimulationData SD;
  SD.readOnlyData  = ROD;
  SD.readWriteData = RWD;
  
  border_print();

  return SD;
}

#define PRNG_SAMPLES_PER_RAY 4
void initialize_ray_kernel(uint64_t base_seed, int ray_id, double length_per_dimension, int n_cells_per_dimension, double inverse_cell_width, RayData RD)
{
    uint64_t offset = ray_id * PRNG_SAMPLES_PER_RAY;
    uint64_t seed = fast_forward_LCG(base_seed, offset);
    RD.location_x[ray_id] = LCG_random_double(&seed) * length_per_dimension;
    RD.location_y[ray_id] = LCG_random_double(&seed) * length_per_dimension;
    //RD.location_z[ray_id] = 0.0;

    // Sample Angle
    double theta = LCG_random_double(&seed) * 2.0 * M_PI;

    // Full 3D Random
    double z = -1.0 + 2.0 * LCG_random_double(&seed);
    if(fabs(z) > 0.99999)
      z = 0.602399;
    
    // TY Quadrature 1 angle
    //double z = 0.602399;
    
    // TY Quadrature 3 angle
    /*
    double polar_sample = LCG_random_double(&seed);
    double TY_polar_quadrature_angles[3] = {0.9860164522440789375086539997532650579039834891670296157543, 0.8431317703366419817667837235644272158155494503259851492417, 0.3599956025898094216886083598423368551278604552721705292320};
    //double TY_polar_quadrature_weights[3] = {0.046233, 0.283619, 0.670148};
    double TY_polar_quadrature_CDF[3] = {0.046233, 0.329852, 1.0};
    int polar_angle;
    for( polar_angle = 0; polar_angle < 3; polar_angle++ )
    {
      if( polar_sample < TY_polar_quadrature_CDF[polar_angle] )
          break;
    }
    double z = TY_polar_quadrature_angles[polar_angle];
    */

    double zo = sqrt(1.0 - z*z);

    // Spherical conversion
    double x = zo * cosf(theta);
    double y = zo * sinf(theta);
    //RD.direction_z[ray_id] = z;

    // Normalize Direction
    double inverse = 1.0 / sqrt( x*x + y*y + z*z );
    x *= inverse;
    y *= inverse;
    z *= inverse;
    
    // Compute Cell ID
    int x_idx = RD.location_x[ray_id] * inverse_cell_width;
    int y_idx = RD.location_y[ray_id] * inverse_cell_width;
    RD.cell_id[ray_id] = y_idx * n_cells_per_dimension + x_idx; 

    RD.direction_x[ray_id] = x;
    RD.direction_y[ray_id] = y;
    
    //if(fabs(z) > 0.999)
    //  printf("ray id: %d direction = [%+.5lf, %+.5lf, %+.5lf] location = [%.5lf, %.5lf] idx = [%d, %d] cell_id = %d\n", ray_id, x, y, z, RD.location_x[ray_id], RD.location_y[ray_id], x_idx, y_idx, RD.cell_id[ray_id]) ;
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
