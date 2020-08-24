#include "minray.h"

void normalize_scalar_flux_kernel(Parameters P, float * new_scalar_flux, int cell, int energy_group)
{
  if( cell >= P.n_cells )
    return;
  if( energy_group >= P.n_energy_groups )
    return;

  new_scalar_flux[(uint64_t) cell * P.n_energy_groups + energy_group] *= P.inverse_total_track_length;
}
