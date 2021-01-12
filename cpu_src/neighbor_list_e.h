#ifndef NEIGHBOR_LIST_E_H
#define NEIGHBOR_LIST_E_H

typedef struct node_{
  int element;
  struct node_ * next;
} Node;

typedef struct{
  Node * head;
  omp_lock_t mutex;
} NeighborList;

typedef struct{
  Node * next;
  int is_finished;
} NeighborListIterator;

void nl_init(         NeighborList * neighborList);
void nl_init_iterator(NeighborList * neighborList, NeighborListIterator * neighborListIterator);
int  nl_read_next(    NeighborList * neighborList, NeighborListIterator * neighborListIterator);
void nl_push_back(    NeighborList * neighborList, int new_elem);

#endif // OPENMC_NEIGHBOR_LIST_A_H
