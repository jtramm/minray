
void nl_push_back(__global int * vectorPool, __global int * vectorPool_idx, int vectorPool_size, __global NeighborList * neighborList, int new_elem)
{
  // Let's first read through the list to make sure it's not there already
  for( int level = 0; level < LEVELS; level++ )
  {
    // Check if this level exists
    int * ptr = atomic_add(&neighborList->ptrs[level],0);
    
    if( ptr == NULL )
    {
      // attempt to push back
    }

  // NOTE: LAST NOTE FOR FRIDAY;
  // I think the scheme will just barely work, but I need to switch to integer indices being stored in ptrs rather than unsigned pointers themselves.
// This is because to avoid duplicate allocations, I will need to change the integers from -1 to -2 first to get permission to allocate, then have the allocate set it to the actual value.
// 
    
    
    
    

  }

  for( int idx = 0; idx < size; idx++ )
  {
    int level = (int) log2(idx+1); 
    int index = idx - (int) pow(2, level-1);

    // this operation needs to be atomic in theory, but I'm not sure you can operate on a ptr in openCL 1.2 atomically. hmmmph.
    // Actually, NO, it's ok! The thread that pushes back will first set ptr to -2, then make allocation, then set pointer to correct value, then set capacity to correct value, then set size to correct value. So, we will not be trying to read from a pointer that is being written to concurrently ever.
    int * ptr = neighborList->ptrs[level];

    // We want to perform an atomic read, but do an add to stay within OpenCL 1.2
    int value = atomic_add(&ptr[index],0);

    if( value == new_elem )
      return;
  }


  // If we have made it here, we will assume our item is not in the list.

  // TODO:
  // Major problem. How to guarantee there are no duplicates?

  // If we FIRST increase the size, this is bad as reads can end up trying to read from the unallocated pointers. Fixable if we do an atomic read there and check for null.
  // If we FIRST increase the size, then we can indeed scan through and guarantee that there are no duplicates. However, this means that capacity, pointers, and the value itself are all not yet written,
  //   - So, we would have to add conditional logic in a bunch of places to check

  // What if we make the size based on -1 content instead of a size variable itself?
  // Then, we just CAS on each element for -1'ness. When we get to the end of a level, we then set the next pointer to -2 and do our thing. Do our thing.

  // I.e., 1) CAS on -2
  // if we are success, we allocate, set all to -1, and then set the pointer to the new allocation.
  // If we are a fail, then just fail to append. This is ok.
  // 
  // One downside to this idea is that if we try to append to a full list, it will get super slow, but who cares about that I guess.

  // Attempt to append the node to the previous one via an atomic compare-and-swap (CAS) operation
  current_node_idx =  atomic_cmpxchg(previous_node_idx, -1, new_node_idx);



}

void nl_init_iterator(__global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  // We want to perform an atomic read, but do an add to stay within OpenCL 1.2
  neighborListIterator->next_idx = atomic_add(&neighborList->head_idx, 0);
}

int nl_read_next(__global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  /*
     typedef struct{
     int * ptrs[8];
     int capacity;
     int size;
     } NeighborList;
   */

  int size;
  size = atomic_add(&neighborList->size, 0);

  int idx = neighborListIterator->idx++;

  if( idx >= size )
    return -1;

  int level = (int) log2(idx+1); 
  int index = idx - (int) pow(2, level-1);

  // this operation needs to be atomic in theory, but I'm not sure you can operate on a ptr in openCL 1.2 atomically. hmmmph.
  // Actually, NO, it's ok! The thread that pushes back will first set ptr to -2, then make allocation, then set pointer to correct value, then set capacity to correct value, then set size to correct value. So, we will not be trying to read from a pointer that is being written to concurrently ever.
  int * ptr = neighborList->ptrs[level];

  // We want to perform an atomic read, but do an add to stay within OpenCL 1.2
  int value = atomic_add(&ptr[index],0);

  return value;
}

// I guess the fundamental problem here is that we can have multiple threads inserting at once. As the 

// Maybe do an atomicCAS on the pointer itself first, and set it to -2. Then, make the allocation knowing that you're the only one that can.


// What about the case where a regular (non-alocating) push back occurs. The push_back moves the size back by 1 atomically, so no other thread will be writing there. However, before it writes, another thread could read, which would be bad. I could make the reads atomic, so they are either -1 or the real value, which would be acceptable.
