#include "minray.h"

void compute_cell_fission_rates_kernel(Parameters P, SimulationData SD, float * scalar_flux, int cell)
{
  if( cell >= P.n_cells )
    return;

  int material_id = SD.readOnlyData.material_id[cell];
  int XS_idx      = material_id * P.n_energy_groups;
  float * nu_Sigma_f = SD.readOnlyData.nu_Sigma_f + XS_idx;

  int flux_idx = cell * P.n_energy_groups;
  scalar_flux += flux_idx;

  float fission_rate = 0.0;
  for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
  {
    fission_rate += nu_Sigma_f[energy_group] * scalar_flux[energy_group];
  }

  SD.readWriteData.cellData.fission_rate[cell] = fission_rate * P.cell_volume;
}
