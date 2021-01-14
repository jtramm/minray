#include<omp.h>
#include<assert.h>
#include"neighbor_list_f.h"

// It is possible for this algorithm to have threads writing to and
// reading from the same location concurrently, resulting in undefined behavior.
// The case where two threads are writing to the list concurrently is protected via
// lock, but the case where one thread is writing and another is reading is not protected.
// The variable in question is the populated length of the statically allocated array.
// As the read from the variable is undefined, any element might be read from. If a bad
// read occurs giving a garbage cell_id, bounds checking can be used on it to ensure
// no segmentation fault occurs. If the bad read produces an in-bounds cell_id, then
// no harm is done other than a wasted CSG check, so the algorithm would still produce
// valid results.

void nl_push_back(NeighborList * neighborList, int new_elem)
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
  assert(idx < NEIGHBOR_SIZE);

  // Write element value to the array
  neighborList->list[idx] = new_elem;

  // NOTE: This write can occur at the same time another thread reads from this variable
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
  omp_init_lock(&neighborList->mutex);
}

int nl_read_next(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  int idx = neighborListIterator->idx++;

  // NOTE: The read of the "length" variable can occur at the same time another thread in nl_push_back is writing to it
  if( idx < neighborList->length && idx < NEIGHBOR_SIZE )
    return neighborList->list[idx];
  else
    return -1;
}
