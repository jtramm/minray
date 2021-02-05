#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_d.h"

// This algorithm can result in a race condition where
// different threads are reading from and writing to the same location
// concurrently in an unprotected manner. The case where
// two threads are appending to the list is protected via lock,
// but the case where one thread is pushing back and another
// is reading from the list is not protected. As the read/write
// conflict can occur over a pointer variable, subsequent use of that
// pointer by the reading thread can result in a seg fault.

void nl_push_back(NeighborListPool neighborListPool, NeighborList * neighborList, int new_elem)
{
  // Lock the object
  omp_set_lock(&neighborList->mutex);

  // If the list is empty (i.e., has no head) add the first item
  if( neighborList->head == NULL )
  {
    neighborList->head = (Node *) malloc(sizeof(Node));
    neighborList->head->element = new_elem;
    neighborList->head->next = NULL;

    // unlock the object and return
    omp_unset_lock(&neighborList->mutex);
    return;
  }

  // If the list is not empty, read through until it is 1) found (do nothing) or 2) the end is reached (append new_elem to it)
  Node * current;
  Node * next = neighborList->head;
  while(next != NULL)
  {
    current = next;
    // If we find the new element is already in the list, we don't want to do anything
    if( current->element == new_elem )
    {
      // unlock the object and return
      omp_unset_lock(&neighborList->mutex);
      return;
    }
    next = current->next;
  }

  // If we have reached the end of this list without finding the object, we need to append it
  current->next = (Node *) malloc(sizeof(Node));
  current->next->element = new_elem;
  current->next->next = NULL;

  // unlock the object
  omp_unset_lock(&neighborList->mutex);
}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->next = neighborList->head;
}

void nl_init(NeighborList * neighborList)
{
  neighborList->head = NULL;
  omp_init_lock(&neighborList->mutex);
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  Node * current = neighborListIterator->next;

  if( current == NULL )
    return -1;

  neighborListIterator->next = current->next;

  return current->element;
}
NeighborListPool nl_init_pool(int n_cells)
{
  NeighborListPool NLP;
  NLP.size = 0;
  NLP.pool = NULL;
  NLP.idx  = NULL;
  return NLP;
}
