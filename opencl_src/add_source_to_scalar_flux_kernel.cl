#include "parameters.h"
__kernel void add_source_to_scalar_flux_kernel(ARGUMENTS)
{
  // Get the index of the current element to be processed
	ulong i = get_global_id(0);
  if( i >= P.n_cells * P.n_energy_groups)
    return;
  ulong cell         = get_group_id(0);
  ulong energy_group = get_local_id(0);

  double Sigma_t_value = Sigma_t[material_id[cell] * P.n_energy_groups + energy_group];

  new_scalar_flux[i] /= (Sigma_t_value * P.cell_volume);
  new_scalar_flux[i] += isotropic_source[i];

  scalar_flux_accumulator[i] += new_scalar_flux[i];
}
