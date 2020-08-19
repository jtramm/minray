#include "minray.h"

int main(int argc, char * argv[])
{
  Parameters P = read_CLI(argc, argv);
  ReadOnlyData ROD = load_2D_C5G7_XS(P);
  SimulationData SD = initialize_simulation(P, ROD);

  return 0;
}
