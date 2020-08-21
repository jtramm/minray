#include "minray.h"

void transport_sweep(Parameters P, SimulationData SD)
{
  // Ray Trace Kernel
  for( int ray = 0; ray < P.n_rays; ray++ )
    ray_trace_kernel(P, SD, SD.readWriteData.rayData, ray);

  // Flux Attenuate Kernel
  for( int ray = 0; ray < P.n_rays; ray++ )
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
      flux_attenuation_kernel(P, SD, ray, energy_group);
}

void print_ray_tracing_buffer(Parameters P, SimulationData SD)
{
  IntersectionData ID = SD.readWriteData.intersectionData;
  for( int r = 0; r < P.n_rays; r++ )
  {
    printf("Ray %d had %d intersections\n", r, ID.n_intersections[r]);
    for( int i = 0; i < ID.n_intersections[r]; i++ )
    {
      int idx = r * P.max_intersections_per_ray + i;
      printf("\tIntersection %d:   cell_id: %d   distance: %.2le   vac reflect: %d\n", i, ID.cell_ids[idx], ID.distances[idx], ID.did_vacuum_reflects[idx]);
    }
  }
}

void update_isotropic_sources(Parameters P, SimulationData SD, double k_eff)
{
  for( int cell = 0; cell < P.n_cells; cell++ )
  {
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
    {
      update_isotropic_sources_kernel(P, SD, cell, energy_group, 1.0/k_eff);
    }
  }
}

void normalize_scalar_flux(Parameters P, SimulationData SD)
{
  for( int cell = 0; cell < P.n_cells; cell++ )
  {
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
    {
      normalize_scalar_flux_kernel(P, SD.readWriteData.cellData.new_scalar_flux, cell, energy_group);
    }
  }
}

void add_source_to_scalar_flux(Parameters P, SimulationData SD)
{
  for( int cell = 0; cell < P.n_cells; cell++ )
  {
    for( int energy_group = 0; energy_group < P.n_energy_groups; energy_group++ )
    {
      add_source_to_scalar_flux_kernel(P, SD, cell, energy_group);
    }
  }
}
  
void compute_cell_fission_rates(Parameters P, SimulationData SD, float * scalar_flux)
{
  for( int cell = 0; cell < P.n_cells; cell++ )
  {
    compute_cell_fission_rates_kernel(P, SD, scalar_flux, cell);
  }
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

double compute_k_eff(Parameters P, SimulationData SD, double old_k_eff)
{
  // Compute old fission rates
  compute_cell_fission_rates(P, SD, SD.readWriteData.cellData.old_scalar_flux);

  // Reduce old fission rates
  double old_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells * P.n_energy_groups);
  
  // Compute new fission rates
  compute_cell_fission_rates(P, SD, SD.readWriteData.cellData.new_scalar_flux);

  // Reduce new fission rates
  double new_total_fission_rate = reduce_sum_float(SD.readWriteData.cellData.fission_rate, P.n_cells * P.n_energy_groups);

  // Update estimate of k-eff
  double new_k_eff = old_k_eff * (new_total_fission_rate / old_total_fission_rate);

  return new_k_eff;
}

void ptr_swap(float ** a, float ** b)
{
  float * tmp = *a;
  *a = *b;
  *b = tmp;
}

void run_simulation(Parameters P, SimulationData SD)
{
  // k is the multiplication factor (or eigenvalue) we are trying to solve for.
  // The eigenvector is the scalar flux vector
  double k_eff = 1.0;

  int active_region = 0;

  for( int iter = 0; iter < P.n_iterations; iter++ )
  {
    // Update Source
    update_isotropic_sources(P, SD, k_eff);

    // Set new scalar fluxes to zero
    memset(SD.readWriteData.cellData.new_scalar_flux, 0, P.n_cells * P.n_energy_groups * sizeof(float));

    // Transport Sweep
    transport_sweep(P, SD);

    // Determine how many FSRs were hit
    int n_cells_hit = reduce_sum_int(SD.readWriteData.cellData.hit_count, P.n_cells);

    // Reset cell hit counters
    memset(SD.readWriteData.cellData.hit_count, 0, P.n_cells * sizeof(int));

    //print_ray_tracing_buffer(P, SD);

    // Normalize Scalar Flux
    normalize_scalar_flux(P, SD);

    // Add Source to Flux
    add_source_to_scalar_flux(P, SD);

    // Compute K-eff
    k_eff = compute_k_eff(P, SD, k_eff);
    
    // Reset scalar flux accumulators if we have finished our inactive iterations
    if( iter >= P.n_inactive_iterations && !active_region )
    {
      active_region = 1;
      memset(SD.readWriteData.cellData.scalar_flux_accumulator, 0, P.n_cells * P.n_energy_groups * sizeof(float));
    }

    // Swap old and new scalar flux pointers
    ptr_swap(&SD.readWriteData.cellData.new_scalar_flux, &SD.readWriteData.cellData.old_scalar_flux);

    // Print status data
    printf("Iter %4d:   k-eff = %.5lf   Percent Cells Missed = %.4lf%%\n", iter, k_eff, (1.0 - (double) n_cells_hit/P.n_cells) * 100.0);
  }
}

int main(int argc, char * argv[])
{
  Parameters P = read_CLI(argc, argv);
  SimulationData SD = initialize_simulation(P);
  initialize_rays(P, SD);
  initialize_fluxes(P, SD);

  run_simulation(P, SD);


  return 0;
}
