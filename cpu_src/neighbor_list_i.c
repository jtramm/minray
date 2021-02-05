#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_i.h"

void nl_push_back(NeighborListPool neighborListPool, NeighborList * neighborList, int new_elem)
{
  // Allocate and initialize the new node for appending
  Node * new_node = (Node *) malloc(sizeof(Node));
  new_node->element = new_elem;
  new_node->next = NULL;

  // We begin with the head node pointer
  Node ** previous_node = &neighborList->head; 
  Node * current_node;

  // Loop to traverse the linked list
  while(1)
  {
    // Attempt to append the node to the previous one via an atomic compare-and-swap (CAS) operation
    // If using OpenMP 5.1, we could just use the included atomic CAS operation, but for now we settle with a builtin.
    __sync_synchronize();
    current_node = __sync_val_compare_and_swap(previous_node, NULL, new_node);
    __sync_synchronize();

    // Result 1: current_node is NULL, meaning that we have reached the end of the list. As such, the CAS resulted in the append working and we are done.
    if( current_node == NULL )
      break;

    // Result 2: current_node is not NULL, so we are not at the end of the list, and the CAS failed to append.
    // If the element of the node the CAS found is the same as what we want to append, we have found a duplicate item so should not try to append.
    // To finish our work, we just need to free the memory we've alloced and then return.
    if( current_node->element == new_elem )
    {
      free(new_node);
      return;
    }

    // Result 3: Similar to result (2), the current_node is not NULL, so we are not at the end of the list and the CAS failed to append.
    // As the element of the node the CAS found is NOT the same as what we want to append, we need to continue on to the next node and try again.
    previous_node = &current_node->next;
  }

}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  #pragma omp atomic read seq_cst
  neighborListIterator->next = neighborList->head;
}

void nl_init(NeighborList * neighborList)
{
  neighborList->head = NULL;
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  Node * current = neighborListIterator->next;

  if( current == NULL )
    return -1;

  #pragma omp atomic read seq_cst
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
