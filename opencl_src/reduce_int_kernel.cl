// This reduction function is adapted in large part from a reduction kernel
// written by Elias Ekondis as part of the cl2-reduce-bench repository available
// online at https://github.com/ekondis/cl2-reduce-bench
// It is distributed under the GPL v2 license.

#define WG_SIZE 256

__kernel void reduce_int_kernel(__global int * in, __global int * result, int size) {
  int gid = get_global_id(0);
  int lid = get_local_id(0);

  // Each thread loads its data into the shared buffer
  __local int localBuffer[WG_SIZE];
  int res = 0;
  if( gid < size)
    res = in[gid];
  localBuffer[lid] = res;
  barrier(CLK_LOCAL_MEM_FENCE);

  // Execute a local memory reduction
  int i = WG_SIZE/2;
  for(; i > 1; i >>= 1)
  {
    if(lid < i)
      localBuffer[lid] = res = res + localBuffer[lid + i];
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  // Perform the wavefront reduction
  for(; i > 0; i >>= 1)
  {
    if(lid < i)
      localBuffer[lid] = res = res + localBuffer[lid + i];
  }

  if(lid == 0 )
    atomic_add(result, res);
}
