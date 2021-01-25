
void nl_push_back(__global Node * nodePool_nodes, __global int * nodePool_idx, int nodePool_size, __global NeighborList * neighborList, int new_elem)
{
  // Pull a fresh node index from the node pool
  int new_node_idx = atomic_inc(nodePool_idx);

  // sanity check (can't assert in OpenCL)
  // assert(new_node_idx < nodePool_size);
  
  // Initialize the new node for appending
  __global Node * new_node = nodePool_nodes + new_node_idx;
  new_node->element = new_elem;
  new_node->next_idx = -1;

  // We begin with the head node pointer
  __global int * previous_node_idx = &neighborList->head_idx; 
  int current_node_idx;

  // Loop to traverse the linked list
  while(1)
  {
    // Attempt to append the node to the previous one via an atomic compare-and-swap (CAS) operation
    current_node_idx =  atomic_cmpxchg(previous_node_idx, -1, new_node_idx);

    // Result 1: current_node is NULL, meaning that we have reached the end of the list. As such, the CAS resulted in the append working and we are done.
    if( current_node_idx == -1 )
      break;

    // Result 2: current_node is not NULL, so we are not at the end of the list, and the CAS failed to append.
    // If the element of the node the CAS found is the same as what we want to append, we have found a duplicate item so should not try to append.
    // To finish our work, we just need to free the memory we've alloced and then return.
    __global Node * current_node = nodePool_nodes + current_node_idx;
    if( current_node->element == new_elem )
    {
      // Note: The current NodePool design allows only for atomic allocation from it, but does not allow for nodes to be given back to it
      // As such, we may end up wasting nodes/memory when duplicates are found.
      //free(new_node);
      return;
    }

    // Result 3: Similar to result (2), the current_node is not NULL, so we are not at the end of the list and the CAS failed to append.
    // As the element of the node the CAS found is NOT the same as what we want to append, we need to continue on to the next node and try again.
    previous_node_idx = &current_node->next_idx;
  }
}

void nl_init_iterator(__global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  // We want to perform an atomic read, but do an add to stay within OpenCL 1.2
  neighborListIterator->next_idx = atomic_add(&neighborList->head_idx, 0);
}

int nl_read_next(__global Node * nodePool_nodes, __global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  __global Node * current = nodePool_nodes + neighborListIterator->next_idx;

  if( neighborListIterator->next_idx == -1 )
    return -1;

  // We want to perform an atomic read, but do an add to stay within OpenCL 1.2
  neighborListIterator->next_idx = atomic_add(&current->next_idx, 0);

  return current->element;
}
