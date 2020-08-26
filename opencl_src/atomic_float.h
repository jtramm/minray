// This floating point atomic function was written by:
// Anca Hamuraru of STREAM HPC on 2/9/2016
// and is available online at:
// https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
void atomicAdd_g_f(volatile __global float *addr, float val)
{
  union {
    unsigned int u32;
    float f32;
  } next, expected, current;
  current.f32 = *addr;
  do {
    expected.f32 = current.f32;
    next.f32 = expected.f32 + val;
    current.u32 = atomic_cmpxchg( (volatile __global unsigned int *)addr,
        expected.u32, next.u32);
  } while( current.u32 != expected.u32 );
}
