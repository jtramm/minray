#include "minray.h"

SimulationResult run_simulation(OpenCLInfo * CL, Parameters P, SimulationData SD)
{
  border_print();
  center_print("SIMULATION", 79);
  border_print();

  double k_eff = 1.0;
  double k_eff_total_accumulator = 0.0;
  double k_eff_sum_of_squares_accumulator = 0.0;

  int is_active_region = 0;

  uint64_t n_total_geometric_intersections = 0;

  double start_time_simulation = get_time();
  double time_in_transport_sweep = 0.0;

  size_t flux_bytes = P.n_cells * P.n_energy_groups * sizeof(float);

  // Power Iteration Loop
  for( int iter = 0; iter < P.n_iterations; iter++ )
  {
    // Reset scalar flux and k-eff accumulators if we have finished our inactive iterations
    if( iter >= P.n_inactive_iterations && !is_active_region )
    {
      is_active_region = 1;
      //printf("resetting scalar flux accumulators...\n");
      clear_array(CL, &SD.readWriteData.cellData.d_scalar_flux_accumulator, flux_bytes);
      //memset(SD.readWriteData.cellData.scalar_flux_accumulator, 0, P.n_cells * P.n_energy_groups * sizeof(float));
      k_eff_total_accumulator = 0.0;
      k_eff_sum_of_squares_accumulator = 0.0;
    }

    // Recompute the isotropic neutron source based on the last iteration's estimate of the scalar flux
    update_isotropic_sources(CL, P, SD, k_eff);

    // Reset this iteration's scalar flux tallies to zero
    //memset(SD.readWriteData.cellData.new_scalar_flux, 0, P.n_cells * P.n_energy_groups * sizeof(float));
    clear_array(CL, &SD.readWriteData.cellData.d_new_scalar_flux, flux_bytes);

    // Run the transport sweep
    double start_time_transport = get_time();
    transport_sweep(CL, P, SD);
    time_in_transport_sweep += get_time() - start_time_transport;

    // Check hit rate to ensure we are running enough rays
    double percent_missed = check_hit_rate(CL, SD, P.n_cells);

    // Normalize the scalar flux tallies to the total distance travelled by all rays this iteration
    normalize_scalar_flux(CL, P, SD);

    // Add the source together with the scalar flux tallies to compute this iteration's estimate of the scalar flux
    add_source_to_scalar_flux(CL, P, SD);

    // Compute a new estimate of the eigenvalue based on the old and new scalar fluxes
    k_eff = compute_k_eff(CL, P, SD, k_eff);
    k_eff_total_accumulator += k_eff;
    k_eff_sum_of_squares_accumulator += k_eff * k_eff;

    // Set old scalar flux to equal the new scalar flux. To optimize, we simply swap the old and new scalar flux pointers
    //ptr_swap(&SD.readWriteData.cellData.new_scalar_flux, &SD.readWriteData.cellData.old_scalar_flux);
    cl_int ret = clEnqueueCopyBuffer(
        CL->command_queue,
        SD.readWriteData.cellData.d_new_scalar_flux,
        SD.readWriteData.cellData.d_old_scalar_flux,
        0,
        0,
        flux_bytes,
        0,
        NULL,
        NULL);
    check(ret);

    // Compute the total number of intersections performed this iteration
    //n_total_geometric_intersections += reduce_sum_int(SD.readWriteData.intersectionData.n_intersections, P.n_rays);
    n_total_geometric_intersections += reduce_intersections(CL, SD, P.n_rays);

    // Output some status data on the results of the power iteration
    print_status_data(iter, k_eff, percent_missed, is_active_region, k_eff_total_accumulator, k_eff_sum_of_squares_accumulator, iter - P.n_inactive_iterations + 1);

  } // End Power Iteration Loop

  double runtime_total = get_time() - start_time_simulation;

  // Gather simulation results
  SimulationResult SR;
  compute_statistics(k_eff_total_accumulator, k_eff_sum_of_squares_accumulator, P.n_active_iterations, &SR.k_eff, &SR.k_eff_std_dev);
  SR.n_geometric_intersections = n_total_geometric_intersections;
  SR.runtime_total = runtime_total;
  SR.runtime_transport_sweep = time_in_transport_sweep;

  // Copy final flux accumulator vector back to the host 
  copy_array_from_device(CL, &SD.readWriteData.cellData.d_scalar_flux_accumulator, SD.readWriteData.cellData.scalar_flux_accumulator, SD.readWriteData.cellData.sz_scalar_flux_accumulator);

  return SR;
}

void update_isotropic_sources(OpenCLInfo * CL, Parameters P, SimulationData SD, double k_eff)
{
  double inv_k_eff = 1.0/k_eff;

  // Update argument
  cl_int ret = clSetKernelArg(CL->kernels.update_isotropic_sources_kernel, 0, sizeof(double), (void *)&inv_k_eff);
  check(ret);

  // Launch kernel
  //printf("Launching update_isotropic_sources kernel...\n");
  size_t global_item_size = P.n_cells * P.n_energy_groups;
  size_t local_item_size = P.n_energy_groups;
  ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.update_isotropic_sources_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);
}

void transport_sweep(OpenCLInfo * CL, Parameters P, SimulationData SD)
{
  size_t global_item_size;
  size_t local_item_size;
  cl_int ret;

  // Launch Ray Tracing kernel
  //printf("Launching ray tracing kernel...\n");
  global_item_size = P.n_rays;
  local_item_size = 8; 
  ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.ray_trace_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);

  // Launch Ray Tracing kernel
  //printf("Launching flux attenuation kernel...\n");
  global_item_size = P.n_rays * P.n_energy_groups;
  local_item_size = P.n_energy_groups;
  ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.flux_attenuation_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);
}


void normalize_scalar_flux(OpenCLInfo * CL, Parameters P, SimulationData SD)
{
  //printf("Launching flux scalar flux normalization kernel...\n");
  size_t global_item_size = P.n_cells * P.n_energy_groups;
  size_t local_item_size = P.n_energy_groups;
  cl_int ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.normalize_scalar_flux_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);
}

void add_source_to_scalar_flux(OpenCLInfo * CL, Parameters P, SimulationData SD)
{
  //printf("Launching add source to scalar flux normalization kernel...\n");
  size_t global_item_size = P.n_cells * P.n_energy_groups;
  size_t local_item_size = P.n_energy_groups;
  cl_int ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.add_source_to_scalar_flux_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);
}

double compute_k_eff(OpenCLInfo * CL, Parameters P, SimulationData SD, double old_k_eff)
{
  // Compute old fission rates
  compute_cell_fission_rates(CL, P, SD, 0.0);

  // Reduce total old fission rate
  //double old_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells * P.n_energy_groups);
  double old_total_fission_rate = reduce_fission_rates(CL, SD, P.n_cells);

  // Compute new fission rates
  compute_cell_fission_rates(CL, P, SD, 1.0);

  // Reduce total new fission rate
  //double new_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells * P.n_energy_groups);
  double new_total_fission_rate = reduce_fission_rates(CL, SD, P.n_cells);

  // Update estimate of k-eff
  double new_k_eff = old_k_eff * (new_total_fission_rate / old_total_fission_rate);

  return new_k_eff;
}

void compute_cell_fission_rates(OpenCLInfo * CL, Parameters P, SimulationData SD, double utility_variable)
{
  cl_int ret = clSetKernelArg(CL->kernels.compute_cell_fission_rates_kernel, 0, sizeof(double), (void *)&utility_variable);
  check(ret);

  //printf("Launching cell fission rates kernel...\n");
  size_t local_item_size = 64;
  size_t global_item_size = ceil(P.n_cells/64.0) * 64.0;
  ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.compute_cell_fission_rates_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);
}

// May need to be a pairwise reduction
double reduce_sum_float(float * a, int size)
{
  double sum = 0.0;

  for( int i = 0; i < size; i++ )
    sum += a[i];

  return sum;
}

float reduce_fission_rates(OpenCLInfo *CL, SimulationData SD, int n_cells)
{
  float sum = 0;
  cl_mem d_sum = copy_array_to_device(CL, CL_MEM_READ_WRITE, (void *) &sum, sizeof(float));

  int argc = 3;
  size_t arg_sz[3];
  void * args[3];

  // Set argument sizes
  arg_sz[0] = sizeof(float *);
  arg_sz[1] = sizeof(float *);
  arg_sz[2] = sizeof(int);

  args[0] = (void *) &SD.readWriteData.cellData.d_fission_rate;
  args[1] = (void *) &d_sum;
  args[2] = (void *) &n_cells;

  set_kernel_arguments(&CL->kernels.reduce_float_kernel, argc, arg_sz, args);

  size_t local_item_size = 64;
  size_t global_item_size = ceil(n_cells/64.0) * 64.0;
  cl_int ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.reduce_float_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);

  copy_array_from_device(CL, &d_sum, (void *) &sum, sizeof(float));

  ret = clReleaseMemObject(d_sum);
  check(ret);

  return sum;
}

int reduce_intersections(OpenCLInfo *CL, SimulationData SD, int n_rays)
{
  int sum = 0;
  cl_mem d_sum = copy_array_to_device(CL, CL_MEM_READ_WRITE, (void *) &sum, sizeof(int));

  int argc = 3;
  size_t arg_sz[3];
  void * args[3];

  // Set argument sizes
  arg_sz[0] = sizeof(int *);
  arg_sz[1] = sizeof(int *);
  arg_sz[2] = sizeof(int);

  args[0] = (void *) &SD.readWriteData.intersectionData.d_n_intersections;
  args[1] = (void *) &d_sum;
  args[2] = (void *) &n_rays;

  set_kernel_arguments(&CL->kernels.reduce_int_kernel, argc, arg_sz, args);

  size_t local_item_size = 64;
  size_t global_item_size = ceil(n_rays/64.0) * 64.0;
  cl_int ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.reduce_int_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);

  copy_array_from_device(CL, &d_sum, (void *) &sum, sizeof(int));

  ret = clReleaseMemObject(d_sum);
  check(ret);

  return sum;
}

int reduce_hit_count(OpenCLInfo *CL, SimulationData SD, int n_cells)
{
  int sum = 0;
  cl_mem d_sum = copy_array_to_device(CL, CL_MEM_READ_WRITE, (void *) &sum, sizeof(int));

  int argc = 3;
  size_t arg_sz[3];
  void * args[3];

  // Set argument sizes
  arg_sz[0] = sizeof(int *);
  arg_sz[1] = sizeof(int *);
  arg_sz[2] = sizeof(int);

  args[0] = (void *) &SD.readWriteData.cellData.d_hit_count;
  args[1] = (void *) &d_sum;
  args[2] = (void *) &n_cells;

  set_kernel_arguments(&CL->kernels.reduce_int_kernel, argc, arg_sz, args);

  size_t local_item_size = 64;
  size_t global_item_size = ceil(n_cells/64.0) * 64.0;
  cl_int ret = clEnqueueNDRangeKernel(CL->command_queue, CL->kernels.reduce_int_kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);

  copy_array_from_device(CL, &d_sum, (void *) &sum, sizeof(int));

  ret = clReleaseMemObject(d_sum);
  check(ret);

  return sum;
}

double check_hit_rate(OpenCLInfo * CL, SimulationData SD, int n_cells)
{
  // Determine how many FSRs were hit
  int n_cells_hit = reduce_hit_count(CL, SD, n_cells);

  // Reset cell hit counters
  //memset(hit_count, 0, n_cells * sizeof(int));
  clear_array(CL, &SD.readWriteData.cellData.d_hit_count, n_cells * sizeof(int));

  // Compute percentage of cells missed
  double percent_missed = (1.0 - (double) n_cells_hit/n_cells) * 100.0;

  return percent_missed;
}
