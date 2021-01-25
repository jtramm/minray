#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include"neighbor_list_j.h"
#include"minray.h"

void nl_init(NeighborList * neighborList)
{
  neighborList->head_idx = -1;
}

NodePool nl_init_nodePool(int n_cells)
{
  Node * nodes = (Node *) malloc(n_cells * AVG_NEIGHBORS_PER_CELL * sizeof(Node));
  int * idx = (int *) malloc(sizeof(int));
  *idx = 0;

  NodePool nodePool;
  nodePool.nodes = nodes;
  nodePool.idx = idx;

  nodePool.sz_nodes = n_cells * AVG_NEIGHBORS_PER_CELL * sizeof(Node);
  nodePool.sz_idx = sizeof(int);
  return nodePool;
}
