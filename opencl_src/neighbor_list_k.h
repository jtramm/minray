#ifndef NEIGHBOR_LIST_K_H
#define NEIGHBOR_LIST_K_H

#define AVG_NEIGHBORS_PER_CELL 11
#define LEVELS 8

typedef struct{
  int ptrs[LEVELS];
} NeighborList;

typedef struct{
  int idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);

#endif // OPENMC_NEIGHBOR_LIST_K_H
