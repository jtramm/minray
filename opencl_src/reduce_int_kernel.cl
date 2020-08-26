// This function was written by Zheming Jin and Hal Finkel
// and is availablle online at:
// https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8425459

/*
#define N 1

__kernel void reduce_int_kernel(__global int * input, global int * output)
{
  __local int sum[N];
  int gid = get_global_id(0);
  int lid = get_local_id(0);
  int WGS = get_local_size(0);

  int16 acc = vload16(gid, input);

  int r = (((acc.s0 + acc.s1) + (acc.s2 + acc.s3)) + 
      ((acc.s4 + acc.s5) + (acc.s6 + acc.s7))) +
    (((acc.s8 + acc.s9) + (acc.sa + acc.sb)) +
     ((acc.sc + acc.sd) + (acc.se + acc.sf)));

  if(lid < N)
    sum[lid] = 0;

  barrier(CLK_LOCAL_MEM_FENCE);
  atomic_add(&sum[gid%N],r);
  barrier(CLK_LOCAL_MEM_FENCE);

  if( lid == WGS - 1 )
    atomic_add(output, sum[0]);
}
*/
//--local_size=1024 --global_size=1024

/*
 * A single group of work items collaborates to
 * perform a sum reduction on the array of floats
 * 'in'.  The reduction result is written to the
 * address 'result'.  The number of elements to
 * be reduced is given by 'size'
 *
 * Written by gpuverify-opencl
 * available online at:
 * https://rise4fun.com/GPUVerify-OpenCL/reduction_correct.cl
 */

/*
#define N 1024

#define tid get_local_id(0)

__kernel void reduce_int_kernel(__global int * in, __global int * result, int size) {

  __local int partial_sums[N];

  partial_sums[tid] = in[tid];
  for(int i = tid + N; i < size; i += N) {
    partial_sums[i] += in[i];
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  for(int d = N/2; d > 0; d >>= 1) {
    if(tid < d) {
      partial_sums[tid] += partial_sums[tid + d];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if(tid == 0) {
    *result = partial_sums[0];
  }
}
*/

// very slow
__kernel void reduce_int_kernel(__global int * in, __global int * result, int size) {
  int gid = get_global_id(0);
  atomic_add(result, in[gid]);
}
