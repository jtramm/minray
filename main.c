#include "minray.h"
  
double get_time(void)
{
  #ifdef _MPI
  return MPI_Wtime();
  #endif

  #ifdef OPENMP
  return omp_get_wtime();
  #endif

  time_t time;
  time = clock();

  return (double) time / (double) CLOCKS_PER_SEC;
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

void print_ray_tracing_buffer(Parameters P, SimulationData SD)
{
  IntersectionData ID = SD.readWriteData.intersectionData;
  for( int r = 0; r < P.n_rays; r++ )
  {
    printf("Ray %d had %d intersections, and is now at location [%.2lf, %.2lf]\n", r, ID.n_intersections[r], SD.readWriteData.rayData.location_x[r], SD.readWriteData.rayData.location_y[r]);
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

void ptr_swap(float ** a, float ** b)
{
  float * tmp = *a;
  *a = *b;
  *b = tmp;
}

void output_thermal_fluxes(Parameters P, SimulationData SD)
{
  char * fname = "thermal_fluxes.dat";
  printf("Writing thermal flux data to file: \"%s\"...\n", fname);
  FILE * fp = fopen(fname, "w");
  int cell_id = 0;
  for( int y = 0; y < P.n_cells_per_dimension; y++)
  {
    for( int x = 0; x < P.n_cells_per_dimension; x++)
    {
      fprintf(fp, "%.3le ", SD.readWriteData.cellData.scalar_flux_accumulator[cell_id++] / P.n_active_iterations);
    }
    fprintf(fp, "\n");
  }
}

void run_simulation(Parameters P, SimulationData SD)
{
  center_print("SIMULATION", 79);
  border_print();
  // k is the multiplication factor (or eigenvalue) we are trying to solve for.
  // The eigenvector is the scalar flux vector
  double k_eff = 1.0;

  int active_region = 0;

  uint64_t n_total_geometric_intersections = 0;

  double start_time = get_time();

  for( int iter = 0; iter < P.n_iterations; iter++ )
  {
    // Update Source
    update_isotropic_sources(P, SD, k_eff);

    // Set new scalar fluxes to zero
    memset(SD.readWriteData.cellData.new_scalar_flux, 0, P.n_cells * P.n_energy_groups * sizeof(float));

    // Transport Sweep
    transport_sweep(P, SD);

    // Check hit rate
    double percent_missed = check_hit_rate(SD.readWriteData.cellData.hit_count, P.n_cells);

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

    // Reduce number of intersections performed
    n_total_geometric_intersections += reduce_sum_int(SD.readWriteData.intersectionData.n_intersections, P.n_rays);

    // Print status data
    if( active_region )
      printf("Iter %5d (Active):     k-eff = %.5lf   Miss Rate = %.4lf%%\n", iter, k_eff, percent_missed);
    else
      printf("Iter %5d (Inactive):   k-eff = %.5lf   Miss Rate = %.4lf%%\n", iter, k_eff, percent_missed);
  }
  
  double end_time = get_time();
  double runtime = end_time - start_time;

  uint64_t n_total_integrations = n_total_geometric_intersections * P.n_energy_groups;
  double time_per_integration = runtime * 1.0e9 / n_total_integrations;

  border_print();
  center_print("RESULTS", 79);
  border_print();

  printf("Time per Integration (TPI) = %.3lf [ns]\n", time_per_integration);
}

int main(int argc, char * argv[])
{
  Parameters P = read_CLI(argc, argv);
  print_user_inputs(P);
  SimulationData SD = initialize_simulation(P);
  initialize_rays(P, SD);
  initialize_fluxes(P, SD);

  run_simulation(P, SD);

  if(P.plotting_enabled)
    plot_3D_vtk(P, SD.readWriteData.cellData.scalar_flux_accumulator, SD.readOnlyData.material_id);

  return 0;
}
