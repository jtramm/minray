#ifndef NEIGHBOR_LIST_B_H
#define NEIGHBOR_LIST_B_H

#define AVG_NEIGHBORS_PER_CELL 11

typedef struct _node NeighborListNode;

typedef struct{
  int list[AVG_NEIGHBORS_PER_CELL];
} NeighborList;

typedef struct{
  int idx;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);
/*
void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator);
int  nl_read_next(    NeighborList * neighborList, NeighborListIterator * neighborListIterator);
void nl_push_back(    NeighborList * neighborList, int new_elem);
*/

#endif // OPENMC_NEIGHBOR_LIST_B_H
