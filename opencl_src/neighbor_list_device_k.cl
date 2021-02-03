
void nl_push_back(__global int * vectorPool, __global int * vectorPool_idx, int vectorPool_size, __global NeighborList * neighborList, int new_elem)
{
//printf("attempting push back of new_elem = %d\n", new_elem);
  // Begin reading through the list
  for( int level = 0; level < LEVELS; level++ )
  {
    // Check if this level exists
    int ptr = atomic_cmpxchg(&neighborList->ptrs[level],-1, -2);
    
    // Case 1 - we read a -1, meaning that this level has not been allocated yet. If this is true, then we we have set the pointer to -2, so we should go ahead and append
    if( ptr == -1 )
    {
      // Push back
      int space = (int) pown(2.0f, level);
      int starting_index = atomic_add(vectorPool_idx, space);
      
      // Set all nodes to -1 (this is done at initialization so we're good)
      
      // Write new element (no other thread can see this so we're protected)
      vectorPool[starting_index] = new_elem;

      // Store new index to list ptr
      atomic_cmpxchg(&neighborList->ptrs[level],-2, starting_index);

  //  printf("push back of %d (w/  alloc) success\n", new_elem);

      return;
    }

    // Case 2 - we read a -2. This means someone else is in the process of appending and allocating, so we should just give up.
    if( ptr == -2 )
{
//printf("push back of %d failed due to lockout\n", new_elem);
      return;
}
    
    // Case 3 - we read a valid index. Thus, we should begin scanning through.
    int space = (int) pown(2.0f, level);  
    for( int i = 0; i < space; i++ )
    {
      // Check to make sure we are not reading off the end of the global array. If so, just stop
      if( ptr + i >= vectorPool_size )
        return;

      // Use a CAS to try to append
      int elem = atomic_cmpxchg(&vectorPool[ptr + i],-1, new_elem);

      // Case A - we read a -1. This means we have successfully inserted the element, so we return
      if( elem == -1 )
      {
//printf("push back of %d (w/o alloc) success\n", new_elem);
        return;
      }

      // Case B - we read a valid value >= 0. This means there is something already there, so we check to see if its a duplicate. If so, we return
      if( elem == new_elem )
{
//printf("push back of %d failed due to duplicate found\n", new_elem);
        return;
}

      // Case C - it is different item, so we simply continue on
    }

    // If we reached this point, then we have a valid pointer, and we scanned all its indices at this level and did not find an opening OR a duplicate, so move onto the next level.
  }
}

/*
    int level = (int) log2(idx+1); 
    int index = idx - (int) pow(2, level-1);
    */

void nl_init_iterator(__global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->idx = 0;
}

int nl_read_next(__global int * vectorPool, int vectorPool_size, __global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  // Get the index we want from our iterator
  int idx = neighborListIterator->idx++;

  // Determine which level and index the iterator corresponds to
  int level = (int) log2((float)(idx+1)); 
  int index = idx - (int) pown(2.0f, level) + 1;
//printf("idx = %d -- level = %d -- level_idx = %d\n", idx, level, index);

  // Atomically read the starting index for this level
  int level_starting_index = atomic_add(&neighborList->ptrs[level], 0);

  // Ensure this level's starting index is actually allocated
  if( level_starting_index < 0 )
    return -1;

  // Determine what the global index to read from the global vector pool is
  int reading_index = level_starting_index + index;

  // Ensure we are not reading off the end of that array for some reason
  if( reading_index >= vectorPool_size )
    return -1;

  // Atomically read from the global vector pool
  int element = atomic_add(&vectorPool[reading_index], 0);

  return element;
}
