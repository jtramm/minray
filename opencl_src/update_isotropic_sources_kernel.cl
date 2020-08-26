#include "parameters.h"

__kernel void update_isotropic_sources_kernel(ARGUMENTS)
{
  // Get the index of the current element to be processed
  ulong i = get_global_id(0);
  if( i >= P.n_cells * P.n_energy_groups)
    return;
  ulong cell            = get_group_id(0);
  ulong energy_group_in = get_local_id(0);

  double inverse_k_eff = utility_variable;

  ulong scalar_flux_idx = cell * P.n_energy_groups;
  ulong XS_base = material_id[cell] * P.n_energy_groups;

  Sigma_s    += XS_base * P.n_energy_groups + energy_group_in * P.n_energy_groups;
  nu_Sigma_f += XS_base;

  old_scalar_flux += scalar_flux_idx;

  float Chi_value =     Chi[    XS_base + energy_group_in];
  float Sigma_t_value = Sigma_t[XS_base + energy_group_in];

  float scatter_source = 0.0;
  float fission_source = 0.0;

  for( int energy_group_out = 0; energy_group_out < P.n_energy_groups; energy_group_out++ )
  {
    scatter_source += Sigma_s[   energy_group_out] * old_scalar_flux[energy_group_out];
    fission_source += nu_Sigma_f[energy_group_out] * old_scalar_flux[energy_group_out];
  }

  fission_source *= Chi_value * inverse_k_eff;
  float new_isotropic_source = (scatter_source + fission_source)  / Sigma_t_value;
  isotropic_source[scalar_flux_idx + energy_group_in] = new_isotropic_source;
}

