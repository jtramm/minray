#include<omp.h>
#include<assert.h>
#include"neighbor_list_b.h"

void nl_push_back(NeighborList * neighborList, int new_elem)
{
  for( int i = 0; i < NEIGHBOR_SIZE; i++)
  {
    // This line checks to see if the new_elem is already in the list
    int retrieved_id;

    // CUDA
    // retrieved_id = atomicCAS(&neighborList->list[i], -1, new_elem);

    // OpenMP 5.1
    #pragma omp atomic compare capture
    {
      retrieved_id = neighborList->list[i];
      if (neighborList->list[i] == -1)
        neighborList->list[i] = new_elem;
    }

    // Case 1: The element was not initialized yet, so the previous line had the effect of setting it to new_elem and returning -1.
    // Case 2: The element was already initialized to the current new_elem, so the atomicCAS call will return new_elem
    if( retrieved_id == -1 || retrieved_id == new_elem)
      return;

    // Case 3: The element was already initialized to a different cell_id, so it will return some other value != -1 and != new_elem
    // so, we continue reading through the list.
  }

  // If we reach here without returning, this means the list is full and NEIGHBOR_SIZE should be increased.
  assert(0);
}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->idx = 0;
}

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < NEIGHBOR_SIZE; i++ )
    neighborList->list[i] = -1; 
}

int nl_read_next(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  int idx = neighborListIterator->idx++;

  if( idx >= NEIGHBOR_SIZE )
    return -1;
  else
  {
    int next_element;
    #pragma omp atomic read
    next_element = neighborList->list[idx];
    return next_element;
  }
}
