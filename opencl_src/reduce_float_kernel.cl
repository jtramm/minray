#include "atomic_float.h"

// This reduction function is adapted in large part from a reduction kernel
// written by Elias Ekondis as part of the cl2-reduce-bench repository available
// online at https://github.com/ekondis/cl2-reduce-bench
// It is distributed under the GPL v2 license available at:
// https://github.com/ekondis/cl2-reduce-bench/blob/master/LICENSE

#define WG_SIZE 256

__kernel void reduce_float_kernel(__global float * in, __global float * result, int size) {
  int gid = get_global_id(0);
  int lid = get_local_id(0);

  // Each thread loads its data into the shared buffer
  __local float localBuffer[WG_SIZE];
  float res = 0;
  if( gid < size)
    res = in[gid];
  localBuffer[lid] = res;
  barrier(CLK_LOCAL_MEM_FENCE);

  // Execute the local memory reduction
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

  if(lid == 0)
    atomicAdd_g_f(result, res);
}
