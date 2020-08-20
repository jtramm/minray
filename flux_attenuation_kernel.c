#include "minray.h"
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

void ray_trace_kernel(Parameters P, SimulationData SD, int ray_id, int energy_group)
{
  float * Q          = GD.D_Q;
  float * phi        = GD.D_new_phi;
  float * ray_fluxes = GD.D_ray_fluxes;
  short * FSR_XS_ids = GD.D_FSR_XS_ids;
  float * Sigma_t    = GD.D_Sigma_t;

  // Each thread within a block represents 1 energy group
  const int e = threadIdx.x;

  // The energy group number
  if( e >= negroups)
    return;

  // The ray ID
  int r = blockIdx.x;
  if( r >= nrays )
    return;

  // The global index that this thred represnts in terms of ray ID and energy group #
  const long r_id = (long) r * (long) negroups + (long) e;

  // The cell ID for the FSR the ray starts life in
  long cell_id = fsr_ids[r*max_intersections];

  // The global index in the cell-wise arrays corresponding to the current FSR's ID and the thread's energy group
  long c_id = cell_id * (long) negroups + e;

  // Set guess for starting angular fluxes to be equal to the FSR's isotropic source
  float ray_flux;
  if( is_fresh )
    ray_flux = Q[c_id];
  else
    ray_flux = ray_fluxes[r_id];

  // Intersection ID offset
  long i_offset = r*max_intersections;

  // Loop over all intersections
  for( int i = 0; i < intersection_count[r]; i++ )
  {
    const long long int intersection_id = i_offset + i;

    // The cell ID for the FSR the ray starts life in
    cell_id = fsr_ids[intersection_id];

    // The global index in the cell-wise arrays corresponding to the current FSR's ID and the thread's energy group
    c_id = cell_id * (long) negroups + e;

    const float dist = distances[intersection_id];

    if( did_vac_reflect[intersection_id] )
      ray_flux = 0.0f;

    int inactive = 1;
    if( i >= active_point[r] )
      inactive = 0;

    // Find XS idx for cell
    const short XS_idx = FSR_XS_ids[cell_id];
    const float sigt = Sigma_t[XS_idx * negroups + e];

    // tau calculation ( tau = Sigma_t * distance )
    const float tau = sigt * dist;

    // Exponential ( exponential = 1 - exp( -tau ) )
    float exponential;
    {
      const float c1n =-1.0000013559236386308f;
      const float c2n =0.23151368626911062025f;
      const float c3n =-0.061481916409314966140f;
      const float c4n =0.0098619906458127653020f;
      const float c5n =-0.0012629460503540849940f;
      const float c6n =0.00010360973791574984608f;
      const float c7n =-0.000013276571933735820960f;

      const float c0d =1.0f;
      const float c1d =-0.73151337729389001396f;
      const float c2d =0.26058381273536471371f;
      const float c3d =-0.059892419041316836940f;
      const float c4d =0.0099070188241094279067f;
      const float c5d =-0.0012623388962473160860f;
      const float c6d =0.00010361277635498731388f;
      const float c7d =-0.000013276569500666698498f;

      float x = -tau;
      float num, den;

      den = c7d;
      den = den * x + c6d;
      den = den * x + c5d;
      den = den * x + c4d;
      den = den * x + c3d;
      den = den * x + c2d;
      den = den * x + c1d;
      den = den * x + c0d;

      num = c7n;
      num = num * x + c6n;
      num = num * x + c5n;
      num = num * x + c4n;
      num = num * x + c3n;
      num = num * x + c2n;
      num = num * x + c1n;
      num = num * x;

      exponential = num / den;
    }

    // Delta_psi = (angular track flux - Q) * exponential
    const float delta_psi =  (ray_flux - Q[c_id]) * exponential;

    if( !inactive )
    {
      // scalar flux += 4 * PI * Delta_psi
      const float tally = FOUR_PI * delta_psi;
      atomicAdd( phi + c_id, tally);
    }

    // track flux -= Delta_Psi
    ray_flux -= delta_psi;

  } // end flux attenuation loop

  // Store fluxes for transmission
  ray_fluxes[r_id] = ray_flux;
}

}
