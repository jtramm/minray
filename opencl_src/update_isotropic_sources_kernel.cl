
__kernel void update_isotropic_sources_kernel(double inverse_k_eff, ulong size, ulong n_energy_groups,
    __global int * material_id,
    __global float * Sigma_s,
    __global float * nu_Sigma_f,
    __global float * scalar_flux,
    __global float * Chi,
    __global float * Sigma_t,
    __global float * isotropic_source,
    __global float * old_scalar_flux
    )
{
  // Get the index of the current element to be processed
  ulong i = get_global_id(0);
  if( i >= size)
    return;
  ulong cell            = get_group_id(0);
  ulong energy_group_in = get_local_id(0);

  ulong scalar_flux_idx = cell * n_energy_groups;
  ulong XS_base = material_id[cell] * n_energy_groups;

  Sigma_s    += XS_base * n_energy_groups + energy_group_in * n_energy_groups;
  nu_Sigma_f += XS_base;

  old_scalar_flux += scalar_flux_idx;

  float Chi_value =     Chi[    XS_base + energy_group_in];
  float Sigma_t_value = Sigma_t[XS_base + energy_group_in];

  float scatter_source = 0.0;
  float fission_source = 0.0;

  for( int energy_group_out = 0; energy_group_out < n_energy_groups; energy_group_out++ )
  {
    scatter_source += Sigma_s[   energy_group_out] * old_scalar_flux[energy_group_out];
    fission_source += nu_Sigma_f[energy_group_out] * old_scalar_flux[energy_group_out];
  }

  fission_source *= Chi_value * inverse_k_eff;
  float new_isotropic_source = (scatter_source + fission_source)  / Sigma_t_value;
  isotropic_source[scalar_flux_idx + energy_group_in] = new_isotropic_source;
}

