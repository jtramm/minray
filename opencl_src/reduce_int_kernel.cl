__kernel void reduce_int_kernel(__global int * in, __global int * result, int size) {
  int gid = get_global_id(0);
  int lid = get_local_id(0);
  int val;
  if(gid >= size)
    val = 0;
  else
    val = in[gid];
  int res = work_group_reduce_add(val);

  if(lid == 0 )
    atomic_add(result, res);
}
