#ifndef NEIGHBOR_LIST_G_H
#define NEIGHBOR_LIST_G_H

#define NEIGHBOR_SIZE 11

typedef struct{
  int list[NEIGHBOR_SIZE];
} NeighborList;

typedef struct{
  int idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);
void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator);
int  nl_read_next(    NeighborList * neighborList, NeighborListIterator * neighborListIterator);
void nl_push_back(    NeighborList * neighborList, int new_elem);

#endif // OPENMC_NEIGHBOR_LIST_B_H
