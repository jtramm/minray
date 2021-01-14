#ifndef NEIGHBOR_LIST_C_H
#define NEIGHBOR_LIST_C_H

typedef struct{
  int * list;
  int length;
  int capacity;
  omp_lock_t mutex;
} NeighborList;

typedef struct{
  int idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);
void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator);
int  nl_read_next(    NeighborList * neighborList, NeighborListIterator * neighborListIterator);
void nl_push_back(    NeighborList * neighborList, int new_elem);

#endif // OPENMC_NEIGHBOR_LIST_C_H
