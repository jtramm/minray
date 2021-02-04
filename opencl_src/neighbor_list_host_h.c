#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include"neighbor_list_h.h"
#include"minray.h"

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < AVG_NEIGHBORS_PER_CELL; i++ )
    neighborList->list[i] = -1; 
}

NeighborListPool nl_pool_init(int n_cells)
{
  NeighborListPool NLP;
  NLP.sz_pool    = 0;
  NLP.sz_idx = 0;
  return NLP;
}
