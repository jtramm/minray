#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include"neighbor_list_k.h"
#include"minray.h"

NeighborListPool nl_pool_init(int n_cells)
{
  NeighborListPool NLP;
  NLP.sz_pool    = n_cells * sizeof(NeighborListNode) * AVG_NEIGHBORS_PER_CELL; 
  NLP.pool       = (NeighborListNode *) malloc(NLP.sz_pool);
  for( int i = 0; i < n_cells * AVG_NEIGHBORS_PER_CELL; i++ )
    NLP.pool[i].element = -1;
  NLP.sz_idx = sizeof(int);
  NLP.idx    = (int *) malloc(NLP.sz_idx);
  *NLP.idx = 0;

  return NLP;
}

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < LEVELS; i++ )
    neighborList->ptrs[i] = -1;
}
