#include<omp.h>
#include<assert.h>
#include"neighbor_list_a.h"

void nl_push_back(NeighborListPool neighborListPool, NeighborList * neighborList, int new_elem)
{
  // Lock the object
  omp_set_lock(&neighborList->mutex);

  // Check to see if the element already exists in the list.
  // If found, return immediately as there's no work to be done.
  for( int i = 0; i < neighborList->length; i++ )
  {
    if( neighborList->list[i] == new_elem )
    {
      omp_unset_lock(&neighborList->mutex);
      return;
    }
  }

  // If we have reached this point, the element is not yet in the list,
  // so we need to add it.

  // Determine the index we want to write to
  int idx = neighborList->length;

  // Ensure we are not writing off the end of the list
  assert(idx < AVG_NEIGHBORS_PER_CELL);

  // Write element value to the array
  neighborList->list[idx] = new_elem;

  // Atomically increase length of list
  #pragma omp atomic seq_cst
  (neighborList->length)++;

  // unlock the object
  omp_unset_lock(&neighborList->mutex);
}

int nl_get_length(NeighborList * neighborList)
{
  int length;

  // Atomically read the length of the list
  #pragma omp atomic read seq_cst
  length = neighborList->length;

  // Ensure we are not reading off the end of the list
  assert(length <= AVG_NEIGHBORS_PER_CELL);

  return length;
}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->idx = 0;
  neighborListIterator->length = nl_get_length(neighborList);
}

void nl_init(NeighborList * neighborList)
{
  neighborList->length = 0;
  omp_init_lock(&neighborList->mutex);
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  int idx = neighborListIterator->idx++;

  if( idx >= neighborListIterator->length )
    return -1;
  else
    return neighborList->list[idx];
}

NeighborListPool nl_init_pool(int n_cells)
{
  NeighborListPool NLP;
  NLP.size = 0;
  NLP.pool = NULL;
  NLP.idx  = NULL;
  return NLP;
}
