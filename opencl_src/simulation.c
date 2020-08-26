#include "minray.h"

void initialize_device_data(Parameters P, SimulationData * SD, OpenCLInfo * CL)
{
  // Copy Read Only Data
  cl_mem_flags mem_type = CL_MEM_READ_ONLY;
  SD->readOnlyData.d_material_id = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.material_id, SD->readOnlyData.sz_material_id);
  SD->readOnlyData.d_nu_Sigma_f  = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.nu_Sigma_f,  SD->readOnlyData.sz_nu_Sigma_f);
  SD->readOnlyData.d_Sigma_f     = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Sigma_f,     SD->readOnlyData.sz_Sigma_f);
  SD->readOnlyData.d_Sigma_t     = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Sigma_t,     SD->readOnlyData.sz_Sigma_t);
  SD->readOnlyData.d_Sigma_s     = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Sigma_s,     SD->readOnlyData.sz_Sigma_s);
  SD->readOnlyData.d_Chi         = copy_array_to_device(CL, mem_type, (void *) SD->readOnlyData.Chi,         SD->readOnlyData.sz_Chi);

  // Copy Read Write Data
  mem_type = CL_MEM_READ_WRITE;
  SD->readWriteData.rayData.angular_flux = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.angular_flux, SD->readWriteData.rayData.sz_angular_flux);
  SD->readWriteData.rayData.location_x   = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.location_x,   SD->readWriteData.rayData.sz_location_x);
  SD->readWriteData.rayData.location_y   = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.location_y,   SD->readWriteData.rayData.sz_location_y);
  SD->readWriteData.rayData.direction_x  = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.direction_x,  SD->readWriteData.rayData.sz_direction_x);
  SD->readWriteData.rayData.direction_y  = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.direction_y,  SD->readWriteData.rayData.sz_direction_y);
  SD->readWriteData.rayData.cell_id      = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.rayData.cell_id,      SD->readWriteData.rayData.sz_cell_id);
  
  SD->readWriteData.intersectionData.n_intersections     = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.n_intersections,     SD->readWriteData.intersectionData.sz_n_intersections);
  SD->readWriteData.intersectionData.cell_ids            = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.cell_ids,            SD->readWriteData.intersectionData.sz_cell_ids);
  SD->readWriteData.intersectionData.distances           = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.distances,           SD->readWriteData.intersectionData.sz_distances);
  SD->readWriteData.intersectionData.did_vacuum_reflects = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.intersectionData.did_vacuum_reflects, SD->readWriteData.intersectionData.sz_did_vacuum_reflects);
  
  SD->readWriteData.cellData.isotropic_source        = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.isotropic_source,        SD->readWriteData.cellData.sz_isotropic_source);
  SD->readWriteData.cellData.new_scalar_flux         = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.new_scalar_flux,         SD->readWriteData.cellData.sz_new_scalar_flux);
  SD->readWriteData.cellData.old_scalar_flux         = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.old_scalar_flux,         SD->readWriteData.cellData.sz_old_scalar_flux);
  SD->readWriteData.cellData.scalar_flux_accumulator = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.scalar_flux_accumulator, SD->readWriteData.cellData.sz_scalar_flux_accumulator);
  SD->readWriteData.cellData.hit_count               = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.hit_count,               SD->readWriteData.cellData.sz_hit_count);
  SD->readWriteData.cellData.fission_rate            = copy_array_to_device(CL, mem_type, (void *) SD->readWriteData.cellData.fission_rate,            SD->readWriteData.cellData.sz_fission_rate);
}

SimulationResult run_simulation(Parameters P, SimulationData SD)
{
  center_print("SIMULATION", 79);
  border_print();

  double k_eff = 1.0;
  double k_eff_total_accumulator = 0.0;
  double k_eff_sum_of_squares_accumulator = 0.0;

  int is_active_region = 0;

  uint64_t n_total_geometric_intersections = 0;

  double start_time_simulation = get_time();
  double time_in_transport_sweep = 0.0;

  // Power Iteration Loop
  for( int iter = 0; iter < P.n_iterations; iter++ )
  {
    // Reset scalar flux and k-eff accumulators if we have finished our inactive iterations
    if( iter >= P.n_inactive_iterations && !is_active_region )
    {
      is_active_region = 1;
      memset(SD.readWriteData.cellData.scalar_flux_accumulator, 0, P.n_cells * P.n_energy_groups * sizeof(float));
      k_eff_total_accumulator = 0.0;
      k_eff_sum_of_squares_accumulator = 0.0;
    }

    // Recompute the isotropic neutron source based on the last iteration's estimate of the scalar flux
    update_isotropic_sources(P, SD, k_eff);

    // Reset this iteration's scalar flux tallies to zero
    memset(SD.readWriteData.cellData.new_scalar_flux, 0, P.n_cells * P.n_energy_groups * sizeof(float));

    // Run the transport sweep
    double start_time_transport = get_time();
    transport_sweep(P, SD);
    time_in_transport_sweep += get_time() - start_time_transport;

    // Check hit rate to ensure we are running enough rays
    double percent_missed = check_hit_rate(SD.readWriteData.cellData.hit_count, P.n_cells);

    // Normalize the scalar flux tallies to the total distance travelled by all rays this iteration
    normalize_scalar_flux(P, SD);

    // Add the source together with the scalar flux tallies to compute this iteration's estimate of the scalar flux
    add_source_to_scalar_flux(P, SD);

    // Compute a new estimate of the eigenvalue based on the old and new scalar fluxes
    k_eff = compute_k_eff(P, SD, k_eff);
    k_eff_total_accumulator += k_eff;
    k_eff_sum_of_squares_accumulator += k_eff * k_eff;

    // Set old scalar flux to equal the new scalar flux. To optimize, we simply swap the old and new scalar flux pointers
    ptr_swap(&SD.readWriteData.cellData.new_scalar_flux, &SD.readWriteData.cellData.old_scalar_flux);

    // Compute the total number of intersections performed this iteration
    n_total_geometric_intersections += reduce_sum_int(SD.readWriteData.intersectionData.n_intersections, P.n_rays);

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

  return SR;
}

void update_isotropic_sources(Parameters P, SimulationData SD, double k_eff)
{
  double inv_k_eff = 1.0/k_eff;

  #pragma omp parallel for
  for( int cell = 0; cell < P.n_cells; cell++ )
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
      update_isotropic_sources_kernel(P, SD, cell, energy_group, inv_k_eff);
}

void transport_sweep(Parameters P, SimulationData SD)
{
  // Ray Trace Kernel
  #pragma omp parallel for
  for( int ray = 0; ray < P.n_rays; ray++ )
    ray_trace_kernel(P, SD, SD.readWriteData.rayData, ray);

  // Flux Attenuate Kernel
  #pragma omp parallel for
  for( int ray = 0; ray < P.n_rays; ray++ )
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
      flux_attenuation_kernel(P, SD, ray, energy_group);
}


void normalize_scalar_flux(Parameters P, SimulationData SD)
{
  #pragma omp parallel for
  for( int cell = 0; cell < P.n_cells; cell++ )
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
      normalize_scalar_flux_kernel(P, SD.readWriteData.cellData.new_scalar_flux, cell, energy_group);
}

void add_source_to_scalar_flux(Parameters P, SimulationData SD)
{
  #pragma omp parallel for
  for( int cell = 0; cell < P.n_cells; cell++ )
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
      add_source_to_scalar_flux_kernel(P, SD, cell, energy_group);
}

double compute_k_eff(Parameters P, SimulationData SD, double old_k_eff)
{
  // Compute old fission rates
  compute_cell_fission_rates(P, SD, SD.readWriteData.cellData.old_scalar_flux);

  // Reduce total old fission rate
  double old_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells * P.n_energy_groups);

  // Compute new fission rates
  compute_cell_fission_rates(P, SD, SD.readWriteData.cellData.new_scalar_flux);

  // Reduce total new fission rate
  double new_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells * P.n_energy_groups);

  // Update estimate of k-eff
  double new_k_eff = old_k_eff * (new_total_fission_rate / old_total_fission_rate);

  return new_k_eff;
}

void compute_cell_fission_rates(Parameters P, SimulationData SD, float * scalar_flux)
{
  #pragma omp parallel for
  for( int cell = 0; cell < P.n_cells; cell++ )
    compute_cell_fission_rates_kernel(P, SD, scalar_flux, cell);
}

// May need to be a pairwise reduction
double reduce_sum_float(float * a, int size)
{
  double sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for( int i = 0; i < size; i++ )
    sum += a[i];

  return sum;
}

int reduce_sum_int(int * a, int size)
{
  int sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for( int i = 0; i < size; i++ )
    sum += a[i];

  return sum;
}


double check_hit_rate(int * hit_count, int n_cells)
{
  // Determine how many FSRs were hit
  int n_cells_hit = reduce_sum_int(hit_count, n_cells);

  // Reset cell hit counters
  memset(hit_count, 0, n_cells * sizeof(int));

  // Compute percentage of cells missed
  double percent_missed = (1.0 - (double) n_cells_hit/n_cells) * 100.0;

  return percent_missed;
}


