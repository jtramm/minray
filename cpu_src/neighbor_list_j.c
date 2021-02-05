#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_j.h"
#include<stdio.h>
#include"atomic_wrapper.h"

void nl_push_back(NeighborListPool neighborListPool, NeighborList * neighborList, int new_elem)
{
  // We begin with the head node pointer
  int * previous_node_idx = &neighborList->head_idx; 
  int current_node_idx;

  // Loop to traverse the linked list
  while(1)
  {
    current_node_idx = atomic_CAS_wrapper(previous_node_idx, -1, -2);

    // Result 1: we read a -1, meaning that we have reached the end of the list. As such, the CAS resulted in updating this element to -2, and we can allocate and append.
    if( current_node_idx == -1 )
    {
      // Pull a fresh node index from the node pool
      int new_node_idx;
      #pragma omp atomic capture seq_cst
      new_node_idx = (*neighborListPool.idx)++;

      assert(new_node_idx < neighborListPool.size);
      
      // Initialize the new node for appending
      NeighborListNode * new_node = neighborListPool.pool + new_node_idx;
      new_node->element = new_elem;
      new_node->next_idx = -1;

      // Write new node onto the list
      #pragma omp atomic write seq_cst
      *previous_node_idx = new_node_idx;

      return;
    }

    // Result 2: current_node is -2, so someone else must be trying to push back and write. We can either
    //    A) Spin on the CAS and try again, or
    //    B) Just give up and fail to append.
    // As we don't care about occasionally failing to append, and want to avoid use of locks, let's just give up and move on.
    if( current_node_idx == -2 )
    {
      return;
    }   

    // Result 3: current_node is a valid index, so we read the node that that index. We find it is storing an element that is the same as what we are trying to append.
    // As such, we just do nothing and return.
    NeighborListNode * current_node = neighborListPool.pool + current_node_idx;
    if( current_node->element == new_elem )
    {
      return;
    }

    // Result 4: If we reach this point, we have found a valid index and checked its element. 
    // As the element of the node the CAS found is NOT the same as what we want to append, we need to continue on to the next node and try again.
    previous_node_idx = &current_node->next_idx;
  }

}

void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  #pragma omp atomic read seq_cst
  neighborListIterator->next_idx = neighborList->head_idx;
}

void nl_init(NeighborList * neighborList)
{
  neighborList->head_idx = -1;
}

int nl_read_next(NeighborListPool neighborListPool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  NeighborListNode * current = neighborListPool.pool + neighborListIterator->next_idx;

  if( neighborListIterator->next_idx < 0 )
    return -1;

  #pragma omp atomic read seq_cst
  neighborListIterator->next_idx = current->next_idx;

  return current->element;
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

