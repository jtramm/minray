#include "parameters.h"

__kernel void normalize_scalar_flux_kernel(ARGUMENTS)
{
  // Get the index of the current element to be processed
	ulong i = get_global_id(0);
  if( i >= P.n_energy_groups * P.n_cells)
    return;

  new_scalar_flux[i] *= P.inverse_total_track_length;
}
