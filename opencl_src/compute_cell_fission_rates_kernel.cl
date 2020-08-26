#include "parameters.h"

__kernel void compute_cell_fission_rates_kernel(ARGUMENTS)
{
  ulong cell = get_global_id(0);
  if( cell >= P.n_cells )
    return;

  __global float * scalar_flux;
  if( utility_variable == 0.0 )
    scalar_flux = old_scalar_flux;
  else
    scalar_flux = new_scalar_flux;
  
  scalar_flux += cell * P.n_energy_groups;
  nu_Sigma_f += material_id[cell] * P.n_energy_groups;

  double fission_rate_accumulator = 0.0;
  for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
  {
    fission_rate_accumulator += nu_Sigma_f[energy_group] * scalar_flux[energy_group];
  }

  fission_rate[cell] = fission_rate_accumulator * P.cell_volume;
}
