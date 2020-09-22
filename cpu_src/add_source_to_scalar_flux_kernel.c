#include "minray.h"

void add_source_to_scalar_flux_kernel(Parameters P, SimulationData SD, int cell, int energy_group)
{
  if( cell >= P.n_cells )
    return;
  if( energy_group >= P.n_energy_groups )
    return;

  int material_id = SD.readOnlyData.material_id[cell];
  RT_FLOAT Sigma_t   = SD.readOnlyData.Sigma_t[material_id * P.n_energy_groups + energy_group];

  float * new_scalar_flux         = SD.readWriteData.cellData.new_scalar_flux;
  float * isotropic_source        = SD.readWriteData.cellData.isotropic_source; 
  float * scalar_flux_accumulator = SD.readWriteData.cellData.scalar_flux_accumulator; 

  uint64_t idx = (uint64_t) cell * P.n_energy_groups + energy_group;

  new_scalar_flux[idx] /= (Sigma_t * P.cell_volume);
  new_scalar_flux[idx] += isotropic_source[idx];

  scalar_flux_accumulator[idx] += new_scalar_flux[idx];
}
