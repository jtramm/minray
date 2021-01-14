#include<omp.h>
#include<assert.h>
#include"neighbor_list_g.h"

void nl_push_back(NeighborList * neighborList, int new_elem)
{
  for( int i = 0; i < NEIGHBOR_SIZE; i++)
  {
    // This line checks to see if the new_elem is already in the list
    int retrieved_id;

    // read value atomically
    #pragma omp atomic read
    retrieved_id = neighborList->list[i];

    // Case 1: If the value is -1, we will then (attempt to) write to it
    // Note: there are no guarantees that the value has not been written to
    // by another thread in the time since we read from it, so we may be overwriting
    // another value or could even be adding a duplicate if another thread overwrote
    // a previous value.
    if(retrieved_id == -1)
    {
      #pragma omp atomic write
      neighborList->list[i] = new_elem;
      return;
    }
    // Case 2: If the value is already equal to the element we want to append, do nothing
    else if(retrieved_id == new_elem)
    {
      return;
    }
    
    // Case 3: The element was already initialized to a different cell_id, so it will return some other value != -1 and != new_elem
    // so, we continue reading through the list.

  }

  // If we reach here without returning, this means the list is full and NEIGHBOR_SIZE should be increased.
  assert(0);
}

// Example of adding the same element twice:
// Thread A reads a -1 and sets that element to 7
// Thread B reads the same -1 and sets the same element to 8
// These both occur at the same time, so they both read -1 before writing. Thread A first sets it to 7 and thread B shortly after sets it to 8.
//
// As the reads and writes are atomic, the value stored by the end must be 7 or 8, so we are weakly correct but have failed to store 7.
//
// Now lets complicate this slightly by saying we have another thread C that is reading through slightly after and wants to store 8.
// It reads this item such that it succesfully reads the 7 that thread A set, so moves onto the next element and stores 8. However,
// Thread B now stores its 8. So, the result is the same cell_ID (8) has now been stored twice. This is still weakly correct.

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
