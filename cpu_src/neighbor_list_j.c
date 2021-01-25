#include<omp.h>
#include<assert.h>
#include<stdlib.h>
#include"neighbor_list_j.h"

void nl_push_back(NodePool * nodePool, NeighborList * neighborList, int new_elem)
{
  // Pull a fresh node index from the node pool
  #pragma omp atomic seq_cst
  int new_node_idx = (*nodePool->idx)++;

  assert(new_node_idx < nodePool->size);
  
  // Initialize the new node for appending
  Node * new_node = nodePool->nodes + new_node_idx;
  new_node->element = new_elem;
  new_node->next = -1;

  // We begin with the head node pointer
  int * previous_node_idx = &neighborList->head_idx; 
  int current_node_idx;

  // Loop to traverse the linked list
  while(1)
  {
    // Attempt to append the node to the previous one via an atomic compare-and-swap (CAS) operation
    // If using OpenMP 5.1, we could just use the included atomic CAS operation, but for now we settle with a builtin.
    __sync_synchronize();
    current_node_idx = __sync_val_compare_and_swap(previous_node_idx, -1, new_node_idx);
    __sync_synchronize();

    // Result 1: current_node is NULL, meaning that we have reached the end of the list. As such, the CAS resulted in the append working and we are done.
    if( current_node_idx == -1 )
      break;

    // Result 2: current_node is not NULL, so we are not at the end of the list, and the CAS failed to append.
    // If the element of the node the CAS found is the same as what we want to append, we have found a duplicate item so should not try to append.
    // To finish our work, we just need to free the memory we've alloced and then return.
    Node * current_node = nodePool->nodes + current_node_idx;
    if( current_node->element == new_elem )
    {
      // Note: The current NodePool design allows only for atomic allocation from it, but does not allow for nodes to be given back to it
      // As such, we may end up wasting nodes/memory when duplicates are found.
      //free(new_node);
      return;
    }

    // Result 3: Similar to result (2), the current_node is not NULL, so we are not at the end of the list and the CAS failed to append.
    // As the element of the node the CAS found is NOT the same as what we want to append, we need to continue on to the next node and try again.
    previous_node_idx = &current_node->next;
  }

}

void nl_init_iterator(NodePool * nodePool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  #pragma omp atomic read seq_cst
  neighborListIterator->next_idx = neighborList->head_idx;
}

void nl_init(Node * nodePool, NeighborList * neighborList)
{
  neighborList->head_idx = -1;
}

int nl_read_next(NodePool * nodePool, NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  Node * current = nodePool->nodes + neighborListIterator->next_idx;

  if( neighborListIterator->next_idx == -1 )
    return -1;

  #pragma omp atomic read seq_cst
  neighborListIterator->next_idx = current->next_idx;

  return current->element;
}

NodePool nl_init_nodePool(int n_cells)
{
  Node * nodes = (Node *) malloc(n_cells * AVG_NEIGHBORS_PER_CELL * sizeof(Node));
  int * idx = (int *) malloc(sizeof(int));
  *idx = 0;

  NodePool nodePool;
  nodePool.nodes = nodes;
  nodePool.idx = idx;
  nodePool.size = n_cells * AVG_NEIGHBORS_PER_CELL;

  return nodePool;
}
