#include "minray.h"

#define FOUR_PI 12.566370614359172953850573533118011536788677597500423283899778369f

void flux_attenuation_kernel(Parameters P, SimulationData SD, int ray_id, int energy_group)
{
  // Cull threads in case of oversubscription
  if( ray_id >= P.n_rays )
    return;
  if( energy_group >= P.n_energy_groups)
    return;

  float * isotropic_source  = SD.readWriteData.cellData.isotropic_source;
  float * new_scalar_flux   = SD.readWriteData.cellData.new_scalar_flux;
  float angular_flux        = SD.readWriteData.rayData.angular_flux[ray_id * P.n_energy_groups + energy_group];

  int * material_id         = SD.readOnlyData.material_id;
  float * Sigma_t           = SD.readOnlyData.Sigma_t;

  int n_intersections       = SD.readWriteData.intersectionData.n_intersections[ray_id];
  int * cell_ids            = SD.readWriteData.intersectionData.cell_ids            + ray_id * P.max_intersections_per_ray;
  double * distances        = SD.readWriteData.intersectionData.distances           + ray_id * P.max_intersections_per_ray;
  int * did_vacuum_reflects = SD.readWriteData.intersectionData.did_vacuum_reflects + ray_id * P.max_intersections_per_ray;

  // Loop over all of this ray's intersections
  for( int i = 0; i < n_intersections; i++ )
  {
    // The cell ID for the FSR the ray starts life in
    int cell_id = cell_ids[i];

    if( did_vacuum_reflects[i] )
      angular_flux = 0.0f;

    // tau calculation ( tau = Sigma_t * distance )
    float tau = Sigma_t[material_id[cell_id] * P.n_energy_groups + energy_group] * distances[i];

    // Exponential ( exponential = 1 - exp( -tau ) )
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

    int flux_idx = cell_id * P.n_energy_groups + energy_group; 

    float delta_psi =  (angular_flux - isotropic_source[flux_idx]) * exponential;

    // TODO NOTE: This operation must be atomic if running in parallel
    #pragma omp atomic
    new_scalar_flux[flux_idx] += FOUR_PI * delta_psi;

    angular_flux -= delta_psi;
    //printf("angular_flux = %.3le\n", angular_flux);

  } // end intersection loop

  // Store final angular flux for next iteration
  SD.readWriteData.rayData.angular_flux[ray_id * P.n_energy_groups + energy_group] = angular_flux;
}
