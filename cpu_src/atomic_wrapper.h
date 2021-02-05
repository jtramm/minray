int atomic_CAS_wrapper(int * ptr, int test, int replace)
{
  // Attempt to append the node to the previous one via an atomic compare-and-swap (CAS) operation
  // If using OpenMP 5.1, we could just use the included atomic CAS operation, but for now we settle with a builtin.
  int val;
  __sync_synchronize();
  val = __sync_val_compare_and_swap(ptr, test, replace);
  __sync_synchronize();
  return val;
}
