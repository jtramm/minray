#include "minray.h"

int main(int argc, char * argv[])
{
  // Read user inputs from command line
  Parameters P = read_CLI(argc, argv);

  // Display inputs and derived inputs
  print_user_inputs(P);

  // Allocate all data required by simulation
  SimulationData SD = initialize_simulation(P);

  // Populate simulation data with starting guesses
  initialize_rays(P, SD);
  initialize_fluxes(P, SD);

  // Run Random Ray Simulation
  SimulationResult SR = run_simulation(P, SD);

  // Display Results
  int is_valid_result = print_results(P, SR);

  // Output VTK plotting file if enabled
  if(P.plotting_enabled)
    plot_3D_vtk(P, SD.readWriteData.cellData.scalar_flux_accumulator, SD.readOnlyData.material_id);

  return is_valid_result;
}
