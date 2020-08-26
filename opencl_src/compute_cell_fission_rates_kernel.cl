__kernel void compute_cell_fission_rates_kernel(__global float * scalar_flux, ulong size, ulong n_energy_groups, double cell_volume,
    __global int * material_id,
    __global float * nu_Sigma_f,
    __global float * fission_rate
    )
{
  ulong cell = get_global_id(0);
  if( cell >= size )
    return;

  nu_Sigma_f += material_id[cell] * n_energy_groups;
  scalar_flux += cell * n_energy_groups;

  double fission_rate_accumulator = 0.0;
  for( int energy_group = 0; energy_group < n_energy_groups; energy_group++ )
  {
    fission_rate_accumulator += nu_Sigma_f[energy_group] * scalar_flux[energy_group];
  }

  fission_rate[cell] = fission_rate_accumulator * cell_volume;
}
