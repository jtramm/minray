#include "minray.h"

SimulationResult run_simulation(Parameters P, SimulationData SD)
{
  center_print("SIMULATION", 79);
  border_print();

  RT_FLOAT k_eff = 1.0;
  RT_FLOAT k_eff_total_accumulator = 0.0;
  RT_FLOAT k_eff_sum_of_squares_accumulator = 0.0;

  int is_active_region = 0;

  uint64_t n_total_geometric_intersections = 0;

  RT_FLOAT start_time_simulation = get_time();
  RT_FLOAT time_in_transport_sweep = 0.0;

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
    RT_FLOAT start_time_transport = get_time();
    transport_sweep(P, SD);
    time_in_transport_sweep += get_time() - start_time_transport;

    // Check hit rate to ensure we are running enough rays
    RT_FLOAT percent_missed = check_hit_rate(SD.readWriteData.cellData.hit_count, P.n_cells);

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
  
  RT_FLOAT runtime_total = get_time() - start_time_simulation;
  
  // Gather simulation results
  SimulationResult SR;
  compute_statistics(k_eff_total_accumulator, k_eff_sum_of_squares_accumulator, P.n_active_iterations, &SR.k_eff, &SR.k_eff_std_dev);
  SR.n_geometric_intersections = n_total_geometric_intersections;
  SR.runtime_total = runtime_total;
  SR.runtime_transport_sweep = time_in_transport_sweep;

  return SR;
}

void update_isotropic_sources(Parameters P, SimulationData SD, RT_FLOAT k_eff)
{
  RT_FLOAT inv_k_eff = 1.0/k_eff;

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

RT_FLOAT compute_k_eff(Parameters P, SimulationData SD, RT_FLOAT old_k_eff)
{
  // Compute old fission rates
  compute_cell_fission_rates(P, SD, SD.readWriteData.cellData.old_scalar_flux);

  // Reduce total old fission rate
  RT_FLOAT old_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells);
  
  // Compute new fission rates
  compute_cell_fission_rates(P, SD, SD.readWriteData.cellData.new_scalar_flux);

  // Reduce total new fission rate
  RT_FLOAT new_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells);

  // Update estimate of k-eff
  RT_FLOAT new_k_eff = old_k_eff * (new_total_fission_rate / old_total_fission_rate);

  return new_k_eff;
}
  
void compute_cell_fission_rates(Parameters P, SimulationData SD, float * scalar_flux)
{
  #pragma omp parallel for
  for( int cell = 0; cell < P.n_cells; cell++ )
    compute_cell_fission_rates_kernel(P, SD, scalar_flux, cell);
}

// May need to be a pairwise reduction
RT_FLOAT reduce_sum_float(float * a, int size)
{
  RT_FLOAT sum = 0.0;

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


RT_FLOAT check_hit_rate(int * hit_count, int n_cells)
{
  // Determine how many FSRs were hit
  int n_cells_hit = reduce_sum_int(hit_count, n_cells);

  // Reset cell hit counters
  memset(hit_count, 0, n_cells * sizeof(int));

  // Compute percentage of cells missed
  RT_FLOAT percent_missed = (1.0 - (RT_FLOAT) n_cells_hit/n_cells) * 100.0;

  return percent_missed;
}


