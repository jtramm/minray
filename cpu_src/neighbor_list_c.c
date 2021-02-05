#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_c.h"

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

  // Resize array dynamically if necessary
  if( idx >= neighborList->capacity )
  {
    if( neighborList->capacity == 0 )
      neighborList->capacity = 1;
    else
      neighborList->capacity *= 2;

    neighborList->list = (int *) realloc(neighborList->list, neighborList->capacity * sizeof(int));
  }

  // Write element value to the array
  neighborList->list[idx] = new_elem;

  // Increase length of list
  (neighborList->length)++;

  // unlock the object
  omp_unset_lock(&neighborList->mutex);
}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->idx = 0;
}

void nl_init(NeighborList * neighborList)
{
  neighborList->length = 0;
  neighborList->capacity = 0;
  neighborList->list = NULL;
  omp_init_lock(&neighborList->mutex);
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  // Lock the object
  omp_set_lock(&neighborList->mutex);

  int idx = neighborListIterator->idx++;

  int element = -1;

  if( idx < neighborList->length )
    element = neighborList->list[idx];
  
  // unlock the object
  omp_unset_lock(&neighborList->mutex);

  return element;
}
NeighborListPool nl_init_pool(int n_cells)
{
  NeighborListPool NLP;
  NLP.size = 0;
  NLP.pool = NULL;
  NLP.idx  = NULL;
  return NLP;
}
