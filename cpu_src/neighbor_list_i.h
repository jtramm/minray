#ifndef NEIGHBOR_LIST_I_H
#define NEIGHBOR_LIST_I_H
typedef struct _neighborlistnode NeighborListNode;

typedef struct{
 NeighborListNode * pool;
 int * idx;
 int size;
} NeighborListPool;

typedef struct node_{
  int element;
  struct node_ * next;
} Node;

typedef struct{
  Node * head;
} NeighborList;

typedef struct{
  Node * next;
} NeighborListIterator;
#endif // OPENMC_NEIGHBOR_LIST_I_H
