#ifndef NEIGHBOR_LIST_K_H
#define NEIGHBOR_LIST_K_H

#define AVG_NEIGHBORS_PER_CELL 11

#ifdef LEVEL_0
#define LEVELS 5
#define FIRST_LEVEL 0
#endif

#ifdef LEVEL_3
#define LEVELS 2
#define FIRST_LEVEL 3
#endif


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

typedef struct{
  NeighborListNode * pool;
  int * idx;
  int size;
} NeighborListPool;


#endif // OPENMC_NEIGHBOR_LIST_K_H
