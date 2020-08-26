__kernel void normalize_scalar_flux_kernel(ulong size, double inverse_total_track_length, __global float * new_scalar_flux)
{
  // Get the index of the current element to be processed
	ulong i = get_global_id(0);
  if( i >= size)
    return;

  new_scalar_flux[i] *= inverse_total_track_length;
}
