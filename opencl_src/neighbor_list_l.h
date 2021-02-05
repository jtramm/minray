#ifndef NEIGHBOR_LIST_L_H
#define NEIGHBOR_LIST_L_H

#define AVG_NEIGHBORS_PER_CELL 11
#define LEVELS 2
#define FIRST_LEVEL 3

typedef struct{
  int element;
} NeighborListNode;

typedef struct{
  int ptrs[LEVELS];
} NeighborList;

typedef struct{
  int level;
  int level_idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);

#endif // OPENMC_NEIGHBOR_LIST_L_H
