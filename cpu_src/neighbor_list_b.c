#include<omp.h>
#include<stdlib.h>
#include<assert.h>
#include"neighbor_list_b.h"

void nl_push_back(NeighborListPool neighborListPool, NeighborList * neighborList, int new_elem)
{
  for( int i = 0; i < AVG_NEIGHBORS_PER_CELL; i++)
  {
    // This line checks to see if the new_elem is already in the list
    int retrieved_id;

    // CUDA
    // retrieved_id = atomicCAS(&neighborList->list[i], -1, new_elem);

    // OpenMP 5.1
    /*
    #pragma omp atomic compare capture
    {
      retrieved_id = neighborList->list[i];
      if (neighborList->list[i] == -1)
        neighborList->list[i] = new_elem;
    }
    */

    // Compiler non-portable builtin atomic CAS
    retrieved_id = __sync_val_compare_and_swap(&neighborList->list[i], -1, new_elem);

    // Case 1: The element was not initialized yet, so the previous line had the effect of setting it to new_elem and returning -1.
    // Case 2: The element was already initialized to the current new_elem, so the atomicCAS call will return new_elem
    if( retrieved_id == -1 || retrieved_id == new_elem)
      return;

    // Case 3: The element was already initialized to a different cell_id, so it will return some other value != -1 and != new_elem
    // so, we continue reading through the list.
  }

  // If we reach here without returning, this means the list is full and AVG_NEIGHBORS_PER_CELL should be increased.
  assert(0);
}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->idx = 0;
}

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < AVG_NEIGHBORS_PER_CELL; i++ )
    neighborList->list[i] = -1; 
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  int idx = neighborListIterator->idx++;

  if( idx >= AVG_NEIGHBORS_PER_CELL )
    return -1;
  else
  {
    int next_element;
    #pragma omp atomic read
    next_element = neighborList->list[idx];
    return next_element;
  }
}
NeighborListPool nl_init_pool(int n_cells)
{
  NeighborListPool NLP;
  NLP.size = 0;
  NLP.pool = NULL;
  NLP.idx  = NULL;
  return NLP;
}
