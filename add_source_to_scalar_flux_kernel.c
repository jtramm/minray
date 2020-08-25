#include "minray.h"

#define FOUR_PI 12.566370614359172953850573533118011536788677597500423283899778369f

void add_source_to_scalar_flux_kernel(Parameters P, SimulationData SD, int cell, int energy_group)
{
  if( cell >= P.n_cells )
    return;
  if( energy_group >= P.n_energy_groups )
    return;

  int material_id = SD.readOnlyData.material_id[cell];
  double Sigma_t   = SD.readOnlyData.Sigma_t[material_id * P.n_energy_groups + energy_group];

  float * new_scalar_flux         = SD.readWriteData.cellData.new_scalar_flux;
  float * isotropic_source        = SD.readWriteData.cellData.isotropic_source; 
  float * scalar_flux_accumulator = SD.readWriteData.cellData.scalar_flux_accumulator; 

  uint64_t idx = (uint64_t) cell * P.n_energy_groups + energy_group;

  new_scalar_flux[idx] /= (Sigma_t * P.cell_volume);
  //new_scalar_flux[idx] /= (Sigma_t);

  //if( !isfinite( new_scalar_flux[idx] ) )
  //  new_scalar_flux[idx] = 0.0f;

  new_scalar_flux[idx] += FOUR_PI * isotropic_source[idx];
  
  if( new_scalar_flux[idx] < 0 )
    new_scalar_flux[idx] = 0;

  //if( new_scalar_flux[idx] > 0.0f )
    scalar_flux_accumulator[idx] += new_scalar_flux[idx];
}
