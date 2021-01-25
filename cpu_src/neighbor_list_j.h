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
  int size;
} NodePool;

typedef struct{
  int head_idx;
} NeighborList;

typedef struct{
  int next_idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);
void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator);
int  nl_read_next(    NodePool * nodePool, NeighborList * neighborList, NeighborListIterator * neighborListIterator);
void nl_push_back(    NodePool * nodePool, NeighborList * neighborList, int new_elem);
NodePool nl_init_nodePool(int n_cells);

#endif // OPENMC_NEIGHBOR_LIST_J_H
