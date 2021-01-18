#include"neighbor_list_b.h"

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < NEIGHBOR_SIZE; i++ )
    neighborList->list[i] = -1; 
}
