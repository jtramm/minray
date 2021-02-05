#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_k.h"
#include<stdio.h>
#include<math.h>

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

void nl_push_back(NeighborListPool neighborListPool, NeighborList * neighborList, int new_elem)
{
  // Begin reading through the list
  for( int level = 0; level < LEVELS; level++ )
  {
    // Check if this level exists
    int ptr = atomic_CAS_wrapper(&neighborList->ptrs[level],-1, -2);

    // Case 1 - we read a -1, meaning that this level has not been allocated yet. If this is true, then we we have set the pointer to -2, so we should go ahead and append
    if( ptr == -1 )
    {
      // Push back
      int space = (int) pow(2.0f, level + FIRST_LEVEL);
      int starting_index;
      #pragma omp atomic capture seq_cst
      starting_index = (*neighborListPool.idx) += space;

      // Set all nodes to -1 (this is done at initialization so we're good)

      // Write new element (no other thread can see this so we're protected)
      neighborListPool.pool[starting_index].element = new_elem;

      // Store new index to list ptr
      atomic_CAS_wrapper(&neighborList->ptrs[level],-2, starting_index);

      return;
    }

    // Case 2 - we read a -2. This means someone else is in the process of appending and allocating, so we should just give up.
    if( ptr == -2 )
    {
      return;
    }

    // Case 3 - we read a valid index. Thus, we should begin scanning through.
    int space = (int) pow(2.0f, level + FIRST_LEVEL);  
    for( int i = 0; i < space; i++ )
    {
      // Check to make sure we are not reading off the end of the global array. If so, just stop
      if( ptr + i >= neighborListPool.size )
        return;

      // Use a CAS to try to append
      int elem = atomic_CAS_wrapper(&neighborListPool.pool[ptr + i].element,-1, new_elem);

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

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->level = 0;
  neighborListIterator->level_idx = 0;
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  // Get the index we want from our iterator
  int level = neighborListIterator->level;
  int idx = neighborListIterator->level_idx++;

  // Check to see if we are reading off the end of the array
  int space = pow(2.0f, level + FIRST_LEVEL);
  if( idx >= space )
  {
    level++;
    neighborListIterator->level++;
    idx = 0;
    neighborListIterator->level_idx = 1;
  }

  // Atomically read the starting index for this level
  int level_starting_index;
  #pragma omp atomic read
  level_starting_index = neighborList->ptrs[level];

  // Ensure this level's starting index is actually allocated
  if( level_starting_index < 0 )
    return -1;

  // Determine what the global index to read from the global vector pool is
  int reading_index = level_starting_index + idx;

  // Ensure we are not reading off the end of that array for some reason
  if( reading_index >= neighborListPool.size )
    return -1;

  // Atomically read from the global vector pool
  int element;
  #pragma omp atomic read seq_cst
  element = neighborListPool.pool[reading_index].element;

  return element;
}

NeighborListPool nl_init_pool(int n_cells)
{
  NeighborListPool NLP;
  NLP.size = n_cells * AVG_NEIGHBORS_PER_CELL;
  size_t sz_pool    = NLP.size * sizeof(NeighborListNode);
  NLP.pool       = (NeighborListNode *) malloc(sz_pool);
  for( int i = 0; i < n_cells * AVG_NEIGHBORS_PER_CELL; i++ )
    NLP.pool[i].element = -1;
  size_t sz_idx = sizeof(int);
  NLP.idx    = (int *) malloc(sz_idx);
  *NLP.idx = 0;

  return NLP;
}

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < LEVELS; i++ )
    neighborList->ptrs[i] = -1;
}
