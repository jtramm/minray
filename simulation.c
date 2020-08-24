#include "minray.h"

SimulationResult run_simulation(Parameters P, SimulationData SD)
{
  center_print("SIMULATION", 79);
  border_print();

  // k is the multiplication factor (or eigenvalue) we are trying to solve for.
  // The eigenvector is the scalar flux vector
  double k_eff = 1.0;

  // We also want to accumulate an average k-eff across all active iterations, and
  // to know the std. dev. of this sample mean.
  double k_eff_total_accumulator = 0.0;
  double k_eff_sum_of_squares_accumulator = 0.0;

  int is_active_region = 0;

  uint64_t n_total_geometric_intersections = 0;

  double start_time = get_time();

  // Power Iteration Loop
  for( int iter = 0; iter < P.n_iterations; iter++ )
  {
    update_isotropic_sources(P, SD, k_eff);

    // Set new scalar fluxes to zero
    memset(SD.readWriteData.cellData.new_scalar_flux, 0, P.n_cells * P.n_energy_groups * sizeof(float));

    transport_sweep(P, SD);

    // Check hit rate to ensure we are running enough rays
    double percent_missed = check_hit_rate(SD.readWriteData.cellData.hit_count, P.n_cells);

    normalize_scalar_flux(P, SD);

    add_source_to_scalar_flux(P, SD);

    k_eff = compute_k_eff(P, SD, k_eff);
    k_eff_total_accumulator += k_eff;
    k_eff_sum_of_squares_accumulator += k_eff * k_eff;

    // Reset scalar flux and k-eff accumulators if we have finished our inactive iterations
    if( iter >= P.n_inactive_iterations && !is_active_region )
    {
      is_active_region = 1;
      memset(SD.readWriteData.cellData.scalar_flux_accumulator, 0, P.n_cells * P.n_energy_groups * sizeof(float));
      k_eff_total_accumulator = 0.0;
      k_eff_sum_of_squares_accumulator = 0.0;
    }

    // Set old scalar flux to equal the new scalar flux. To optimize, we simply swap the old and new scalar flux pointers
    ptr_swap(&SD.readWriteData.cellData.new_scalar_flux, &SD.readWriteData.cellData.old_scalar_flux);

    // Compute the total number of intersections performed this iteration
    n_total_geometric_intersections += reduce_sum_int(SD.readWriteData.intersectionData.n_intersections, P.n_rays);

    print_status_data(iter, k_eff, percent_missed, is_active_region);

  } // End Power Iteration Loop
  
  double runtime = get_time() - start_time;
  
  // Gather simulation results
  SimulationResult SR;
  int n = P.n_active_iterations;
  SR.k_eff_std_dev = sqrt( (k_eff_sum_of_squares_accumulator - k_eff_total_accumulator * k_eff_total_accumulator / n ) / n);
  SR.k_eff_std_dev /= sqrt(n);
  SR.k_eff = k_eff_total_accumulator / P.n_active_iterations;

  SR.n_geometric_intersections = n_total_geometric_intersections;
  SR.runtime = runtime;

  return SR;
}

void update_isotropic_sources(Parameters P, SimulationData SD, double k_eff)
{
  double inv_k_eff = 1.0/k_eff;

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
  for( int cell = 0; cell < P.n_cells; cell++ )
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
      normalize_scalar_flux_kernel(P, SD.readWriteData.cellData.new_scalar_flux, cell, energy_group);
}

void add_source_to_scalar_flux(Parameters P, SimulationData SD)
{
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
  for( int cell = 0; cell < P.n_cells; cell++ )
    compute_cell_fission_rates_kernel(P, SD, scalar_flux, cell);
}

double reduce_sum_float(float * a, int size)
{
  double sum = 0.0;

  for( int i = 0; i < size; i++ )
    sum += a[i];

  return sum;
}

int reduce_sum_int(int * a, int size)
{
  int sum = 0.0;

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

