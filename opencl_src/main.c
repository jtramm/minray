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
  
  center_print("DEVICE INITIALIZATION", 79);
  border_print();
  // Initialize OpenCL Device
  OpenCLInfo CL = initialize_device();

  // Allocate and move data to device
  initialize_device_data(&SD, &CL);

  // Load and compile OpenCL kernels
  initialize_kernels(&CL);

  // Load static kernel arguments into kernel
  load_kernel_arguments(&P, &SD, &CL);

  // Run Random Ray Simulation
  SimulationResult SR = run_simulation(&CL, P, SD);

  // Display Results
  int is_valid_result = print_results(P, SR);

  // Output VTK plotting file if enabled
  if(P.plotting_enabled)
    plot_3D_vtk(P, SD.readWriteData.cellData.scalar_flux_accumulator, SD.readOnlyData.material_id);

  return is_valid_result;
}
