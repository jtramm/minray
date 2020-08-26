#include "minray.h"

void update_isotropic_sources_kernel(Parameters P, SimulationData SD, int cell, int energy_group_in, double inverse_k_eff)
{
  // Cull threads if oversubscribed
  if( cell >= P.n_cells )
    return;
  if( energy_group_in >= P.n_energy_groups )
    return;

  int material_id = SD.readOnlyData.material_id[cell];

  const uint64_t scalar_flux_idx = cell * P.n_energy_groups;
  const int XS_base = material_id * P.n_energy_groups;

  const float * Sigma_s = SD.readOnlyData.Sigma_s + XS_base * P.n_energy_groups + energy_group_in * P.n_energy_groups;
  const float * nu_Sigma_f = SD.readOnlyData.nu_Sigma_f + XS_base;

  const float * scalar_flux = SD.readWriteData.cellData.old_scalar_flux + scalar_flux_idx;
  
  float Chi =     SD.readOnlyData.Chi[    XS_base + energy_group_in];
  float Sigma_t = SD.readOnlyData.Sigma_t[XS_base + energy_group_in];

  float scatter_source = 0.0;
  float fission_source = 0.0;

  for( int energy_group_out = 0; energy_group_out < P.n_energy_groups; energy_group_out++ )
  {
    scatter_source += Sigma_s[   energy_group_out] * scalar_flux[energy_group_out];
    fission_source += nu_Sigma_f[energy_group_out] * scalar_flux[energy_group_out];
  }

  fission_source *= Chi * inverse_k_eff;
  float new_isotropic_source = (scatter_source + fission_source)  / Sigma_t;
  SD.readWriteData.cellData.isotropic_source[scalar_flux_idx + energy_group_in] = new_isotropic_source;
}

