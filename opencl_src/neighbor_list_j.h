#ifndef NEIGHBOR_LIST_J_H
#define NEIGHBOR_LIST_J_H

#define AVG_NEIGHBORS_PER_CELL 11

typedef struct node_{
  int element;
  int next_idx;
} Node;

typedef struct{
  Node * nodes;
  int * idx;
  cl_mem d_nodes;
  size_t sz_nodes;
  cl_mem d_idx;
  size_t sz_idx;
} NodePool;

typedef struct{
  int element;
  int next_idx;
} NeighborList;

typedef struct{
  int next_idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);
NodePool nl_init_nodePool(int n_cells);

#endif // OPENMC_NEIGHBOR_LIST_J_H
