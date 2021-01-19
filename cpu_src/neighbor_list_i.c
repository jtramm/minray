#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_i.h"

void nl_push_back(NeighborList * neighborList, int new_elem)
{
  // Allocate and initialize the new node for appending
  Node * new_node = (Node *) malloc(sizeof(Node));
  new_node->element = new_elem;
  new_node->next = NULL;

  // Attempt to push back the new node to the list.
  //Node * current_node = atomicCAS(neighborList->head, NULL, new_node);
  Node * current_node;
  __sync_synchronize();
  current_node = __sync_val_compare_and_swap(&neighborList->head, NULL, new_node);
  __sync_synchronize();

  // Result 1: current_node is NULL, meaning that no elements were present in the list, so our append succeeded and we are all done.
  if( current_node == NULL )
    return;

  // Result 2: current_node is not NULL, so the append failed. If the element is the same as what we want to append, we just need to free the memory we've alloced and then are done
  if( current_node->element == new_elem )
  {
    free(new_node);
    return;
  }

  // Result 3: current_node is not NULL, so the append failed. If the element is not what we wanted to append, we need to continue on to the next node
  // if( current_node->element != new_elem ) // (implicit)

  while(1)
  {
    __sync_synchronize();
    current_node = __sync_val_compare_and_swap(&current_node->next, NULL, new_node);
    __sync_synchronize();

    // Result 1: current_node is NULL, meaning that no elements were present in the list, so our append succeeded and we are all done.
    if( current_node == NULL )
      break;

    // Result 2: current_node is not NULL, so the append failed. If the element is the same as what we want to append, we just need to free the memory we've alloced and then are done
    if( current_node->element == new_elem )
    {
      free(new_node);
      return;
    }
    // Result 3: current_node is not NULL, so the append failed. If the element is not what we wanted to append, we need to continue on to the next node
  }

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

int nl_read_next(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  Node * current = neighborListIterator->next;

  if( current == NULL )
    return -1;

  #pragma omp atomic read seq_cst
  neighborListIterator->next = current->next;

  return current->element;
}
