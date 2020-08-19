#include "minray.h"

int main(int argc, char * argv[])
{
  Parameters P = read_CLI(argc, argv);
  ReadOnlyData R = load_2D_C5G7_XS(P);

  return 0;
}
