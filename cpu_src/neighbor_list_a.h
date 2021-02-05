#ifndef NEIGHBOR_LIST_A_H
#define NEIGHBOR_LIST_A_H

#define AVG_NEIGHBORS_PER_CELL 11

typedef struct _neighborlistnode NeighborListNode;

typedef struct{
 NeighborListNode * pool;
 int * idx;
 int size;
} NeighborListPool;

typedef struct{
  int list[AVG_NEIGHBORS_PER_CELL];
  int length;
  omp_lock_t mutex;
} NeighborList;

typedef struct{
  int idx;
  int length;
} NeighborListIterator;

#endif // OPENMC_NEIGHBOR_LIST_A_H
