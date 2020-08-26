#include "atomic_float.h"

// very slow
__kernel void reduce_float_kernel(__global float * in, __global float * result, int size) {
  int gid = get_global_id(0);
  atomicAdd_g_f(result, in[gid]);
}
