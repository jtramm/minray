#ifndef NEIGHBOR_LIST_C_H
#define NEIGHBOR_LIST_C_H
typedef struct _neighborlistnode NeighborListNode;

typedef struct{
 NeighborListNode * pool;
 int * idx;
 int size;
} NeighborListPool;

typedef struct{
  int * list;
  int length;
  int capacity;
  omp_lock_t mutex;
} NeighborList;

typedef struct{
  int idx;
} NeighborListIterator;

#endif // OPENMC_NEIGHBOR_LIST_C_H
