#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include"neighbor_list_k.h"
#include"minray.h"

void nl_init(NeighborList * neighborList)
{
  for( int i = 0; i < LEVELS; i++ )
    neighborList->ptrs[i] = -1;
}
