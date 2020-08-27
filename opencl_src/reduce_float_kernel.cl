#include "atomic_float.h"

// very slow
__kernel void reduce_float_kernel(__global float * in, __global float * result, int size) {
  int gid = get_global_id(0);
  int lid = get_local_id(0);
float val;
if(gid >= size)
val = 0;
else
  val = in[gid];
  float res = work_group_reduce_add(val);

  if(lid == 0 )
    atomicAdd_g_f(result, res);
  //atomic_add(result, in[gid]);
}
