#include "parameters.h"
#include "atomic_float.h"

__kernel void flux_attenuation_kernel(ARGUMENTS)
{
  // Get the index of the current element to be processed
  ulong ray_group_idx = get_global_id(0);
  if( ray_group_idx >= P.n_rays * P.n_energy_groups)
    return;
  ulong ray_id       = get_group_id(0);
  ulong energy_group = get_local_id(0);

  // Indexing
  float ray_angular_flux    = angular_flux[ray_group_idx];
  int ray_n_intersections   = n_intersections[ray_id];

  ulong ray_offset = ray_id * P.max_intersections_per_ray;
  cell_ids            += ray_offset;
  distances           += ray_offset;
  did_vacuum_reflects += ray_offset;

  // Loop over all of this ray's intersections
  for( int i = 0; i < ray_n_intersections; i++ )
  {
    // The cell ID for the FSR the ray starts life in
    ulong cell_id = cell_ids[i];

    if( did_vacuum_reflects[i] )
      ray_angular_flux = 0.0f;

    // tau calculation ( tau = Sigma_t * distance )
    float tau = Sigma_t[material_id[cell_id] * P.n_energy_groups + energy_group] * distances[i];

    /////////////////////////////////////////////////////////////////////
    // Exponential Computation ( exponential = 1 - exp( -tau ) )
    /////////////////////////////////////////////////////////////////////

    // Intrinsic version:
    // float exponential = -expm1(-tau);

    // Explicit version:
    float exponential;
    {
      const float c1n =-1.0000013559236386308f;
      const float c2n = 0.23151368626911062025f;
      const float c3n =-0.061481916409314966140f;
      const float c4n = 0.0098619906458127653020f;
      const float c5n =-0.0012629460503540849940f;
      const float c6n = 0.00010360973791574984608f;
      const float c7n =-0.000013276571933735820960f;

      const float c0d = 1.0f;
      const float c1d =-0.73151337729389001396f;
      const float c2d = 0.26058381273536471371f;
      const float c3d =-0.059892419041316836940f;
      const float c4d = 0.0099070188241094279067f;
      const float c5d =-0.0012623388962473160860f;
      const float c6d = 0.00010361277635498731388f;
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
    /////////////////////////////////////////////////////////////////////

    ulong flux_idx = cell_id * P.n_energy_groups + energy_group; 

    float delta_psi = (ray_angular_flux - isotropic_source[flux_idx]) * exponential;

    //#pragma omp atomic
    //new_scalar_flux[flux_idx] += delta_psi;
    atomicAdd_g_f(new_scalar_flux + flux_idx, delta_psi);

    ray_angular_flux -= delta_psi;

  } // end intersection loop

  // Store final angular flux for next iteration
  angular_flux[ray_group_idx] = ray_angular_flux;
}
