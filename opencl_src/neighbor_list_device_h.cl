
void nl_push_back(__global NeighborList * neighborList, int new_elem)
{
  for( int i = 0; i < NEIGHBOR_SIZE; i++)
  {
    // This line checks to see if the new_elem is already in the list
    int retrieved_id = neighborList->list[i];

    // Case 1: The element was not initialized yet, so the previous line had the effect of setting it to new_elem and returning -1.
    if(retrieved_id == -1)
    {
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
  //assert(0);
}

void nl_init_iterator(__global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  neighborListIterator->idx = 0;
}

int nl_read_next(__global NeighborList * neighborList, NeighborListIterator * neighborListIterator)
{
  int idx = neighborListIterator->idx++;

  if( idx >= NEIGHBOR_SIZE )
    return -1;
  else
  {
    int next_element = neighborList->list[idx];
    return next_element;
  }
}
