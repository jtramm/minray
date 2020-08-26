__kernel void add_source_to_scalar_flux_kernel(ulong size, ulong n_energy_groups, double cell_volume,
    __global int * material_id,
    __global double * Sigma_t,
    __global float * new_scalar_flux,
    __global float * isotropic_source,
    __global float * scalar_flux_accumulator)
{
  // Get the index of the current element to be processed
	ulong i = get_global_id(0);
  if( i >= size)
    return;
  ulong cell         = get_group_id(0);
  ulong energy_group = get_local_id(0);

  double Sigma_t_value = Sigma_t[material_id[cell] * n_energy_groups + energy_group];

  new_scalar_flux[i] /= (Sigma_t_value * cell_volume);
  new_scalar_flux[i] += isotropic_source[i];

  scalar_flux_accumulator[i] += new_scalar_flux[i];
}
