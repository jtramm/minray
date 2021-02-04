#ifndef NEIGHBOR_LIST_J_H
#define NEIGHBOR_LIST_J_H

#define AVG_NEIGHBORS_PER_CELL 11

typedef struct{
  int element;
  int next_idx;
} NeighborListNode;

typedef struct{
  int head_idx;
} NeighborList;

typedef struct{
  int next_idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);

#endif // OPENMC_NEIGHBOR_LIST_J_H
