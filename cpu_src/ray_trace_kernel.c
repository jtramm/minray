#include "minray.h"

void cartesian_ray_trace_pol(RT_FLOAT* x_actual, RT_FLOAT* y_actual, RT_FLOAT cell_width, int x_idx, int y_idx, RT_FLOAT x_dir, RT_FLOAT y_dir);

void ray_trace_kernel(Parameters P, SimulationData SD, RayData rayData, uint64_t ray_id)
{
  RT_FLOAT distance_travelled = 0.0;
  int intersection_id = 0;

  RT_FLOAT x =       rayData.location_x[  ray_id];
  RT_FLOAT y =       rayData.location_y[  ray_id];
  RT_FLOAT x_dir =   rayData.direction_x[ ray_id];
  RT_FLOAT y_dir =   rayData.direction_y[ ray_id];
  int cell_id =      rayData.cell_id[     ray_id];
  int last_surface = rayData.last_surface[ray_id];

  int x_idx = cell_id % P.n_cells_per_dimension;
  int y_idx = cell_id / P.n_cells_per_dimension;

  int just_hit_vacuum = 0;
  int is_terminal = 0;

  // We run this loop until either:
  // 1) The maximum number of intersections has been reached (not typical -- would indicate an error)
  // 2) The ray has reached its set distance (typical operation)
  for( intersection_id = 0; (intersection_id < P.max_intersections_per_ray) && (distance_travelled < P.distance_per_ray); intersection_id++ )
  {
    // Perform ray trace through a Cartesian geometry
    TraceResult trace = cartesian_ray_trace(x, y, P.cell_width, x_idx, y_idx, x_dir, y_dir, last_surface);
  
    // Check to see if ray has reached its maximum distance. Truncate if needed
    if(distance_travelled + trace.distance_to_surface >= P.distance_per_ray)
    {
      trace.distance_to_surface = (P.distance_per_ray - distance_travelled) + BUMP;
      is_terminal = 1;
    }
  
    // Record intersection information for use by flux attenuation kernel
    uint64_t global_intersection_id = ray_id * P.max_intersections_per_ray + intersection_id;
    SD.readWriteData.intersectionData.distances[          global_intersection_id] = trace.distance_to_surface;
    SD.readWriteData.intersectionData.cell_ids[           global_intersection_id] = cell_id;
    SD.readWriteData.intersectionData.did_vacuum_reflects[global_intersection_id] = just_hit_vacuum;
    SD.readWriteData.cellData.hit_count[                                 cell_id] = 1;
    just_hit_vacuum = 0;

    // Move ray forward to intersection surface
    x += x_dir * trace.distance_to_surface;
    y += y_dir * trace.distance_to_surface;
    
    // Create a test point inside the next cell
    RT_FLOAT x_across_surface = x + trace.surface_normal_x * BUMP;
    RT_FLOAT y_across_surface = y + trace.surface_normal_y * BUMP;
    
    // Look up the "neighbor" cell id of the test point.
    // This function also gives us some info on if we hit a boundary, and what type it was.
    //CellLookup lookup = find_cell_id(P, x_across_surface, y_across_surface);
    CellLookup lookup;
    if( is_terminal )
    {
      lookup.cell_id = cell_id;
      lookup.cartesian_cell_idx_x = x_idx;
      lookup.cartesian_cell_idx_y = y_idx;
      lookup.boundary_condition = NONE;
    }
    else
      lookup = find_cell_id(P, x_across_surface, y_across_surface, trace.last_surface, x_idx, y_idx);

    if(!(lookup.cell_id != cell_id || is_terminal) )
    {
      //printf("distance = %f BC = %d current cell_id = %d neighbor id = %d\n", trace.distance_to_surface, lookup.boundary_condition, cell_id, lookup.cell_id);
      printf("distance = %f BC = %d x_norm = %f y_norm = %f\n", trace.distance_to_surface, lookup.boundary_condition, trace.surface_normal_x, trace.surface_normal_y);
      print_ray(x, y, x_dir, y_dir, cell_id);
    }
    
    // A sanity check
    assert(lookup.cell_id != cell_id || is_terminal);

    // If we hit an outer boundary, reflect the ray
    if( lookup.boundary_condition != NONE && !is_terminal )
    {
      trace.surface_normal_x *= -1.0;
      trace.surface_normal_y *= -1.0;
      if( trace.surface_normal_x )
        x_dir *= -1.0;
      else
        y_dir *= -1.0;
      trace.last_surface *= -1;
    }

    // Note if we hit a vacuum boundary
    if( lookup.boundary_condition == VACUUM )
      just_hit_vacuum = 1;

    // If we didn't hit a boundary, the ray is moved into the next cell
    if( lookup.boundary_condition == NONE )
    {
      cell_id = lookup.cell_id;
      x_idx =   lookup.cartesian_cell_idx_x;
      y_idx =   lookup.cartesian_cell_idx_y;
    }

    // Move ray off of surface
    /*
    x += trace.surface_normal_x * BUMP;
    y += trace.surface_normal_y * BUMP;
    */

    last_surface = trace.last_surface;

    // Add this intersection's distance to the total for the ray
    distance_travelled += trace.distance_to_surface;

    // Some sanity checks (can be disabled if desired)
    assert(cell_id >= 0 && cell_id < P.n_cells);
    //printf("distance = %f BC = %d\n", trace.distance_to_surface, lookup.boundary_condition);
    //print_ray(x, y, x_dir, y_dir, cell_id);
    /*
    if( !(x > 0.0 && y > 0.0 && x < P.length_per_dimension && y < P.length_per_dimension))
    {
      printf("distance = %f BC = %d x_norm = %f y_norm = %f\n", trace.distance_to_surface, lookup.boundary_condition, trace.surface_normal_x, trace.surface_normal_y);
      print_ray(x, y, x_dir, y_dir, cell_id);
    }

    assert(x > 0.0 && y > 0.0 && x < P.length_per_dimension && y < P.length_per_dimension);
    */
  }

  if(intersection_id >= P.max_intersections_per_ray)
  {
    printf("WARNING: Increase max number of intersections per ray\n");
    print_ray(x, y, x_dir, y_dir, cell_id);
  }
  
  // Bank the ray's status for use in the next iteration
  rayData.location_x[ ray_id] = x;
  rayData.location_y[ ray_id] = y;
  rayData.direction_x[ray_id] = x_dir;
  rayData.direction_y[ray_id] = y_dir;
  rayData.cell_id[    ray_id] = cell_id;
  //rayData.last_surface[    ray_id] = last_surface;
  rayData.last_surface[    ray_id] = 0;
    
  // Bank number of intersections that this ray had this iteration
  SD.readWriteData.intersectionData.n_intersections[ray_id] = intersection_id;
}

CellLookup find_cell_id(Parameters P, RT_FLOAT x, RT_FLOAT y, int last_surface, int idx_x, int idx_y)
{
  int cartesian_cell_idx_x;
  int cartesian_cell_idx_y;

  if( last_surface == 0 )
  {
    cartesian_cell_idx_x = floor(x * P.inverse_cell_width);
    cartesian_cell_idx_y = floor(y * P.inverse_cell_width);
  }
  else if( last_surface == 1 )
  {
    cartesian_cell_idx_x = idx_x + 1;
    cartesian_cell_idx_y = idx_y;
  }
  else if( last_surface == 2 )
  {
    cartesian_cell_idx_x = idx_x;
    cartesian_cell_idx_y = idx_y + 1;
  }
  else if( last_surface == -1 )
  {
    cartesian_cell_idx_x = idx_x - 1;
    cartesian_cell_idx_y = idx_y;
  }
  else if( last_surface == -2 )
  {
    cartesian_cell_idx_x = idx_x;
    cartesian_cell_idx_y = idx_y - 1;
  }
  int boundary_x = 1;
  int boundary_y = 1;

  if( cartesian_cell_idx_x >= P.n_cells_per_dimension )
    boundary_x = 2;
  else if( cartesian_cell_idx_x < 0 )
    boundary_x = 0;
  else if( cartesian_cell_idx_y >= P.n_cells_per_dimension )
    boundary_y = 2;
  else if( cartesian_cell_idx_y < 0 )
    boundary_y = 0;

  int boundary_condition = P.boundary_conditions[boundary_x][boundary_y];

  int cell_id = cartesian_cell_idx_y * P.n_cells_per_dimension + cartesian_cell_idx_x; 

  CellLookup lookup;
  lookup.cell_id = cell_id;
  lookup.cartesian_cell_idx_x = cartesian_cell_idx_x;
  lookup.cartesian_cell_idx_y = cartesian_cell_idx_y;
  lookup.boundary_condition = boundary_condition;
  return lookup;
}

TraceResult cartesian_ray_trace(RT_FLOAT x, RT_FLOAT y, RT_FLOAT cell_width, int x_idx, int y_idx, RT_FLOAT x_dir, RT_FLOAT y_dir, int last_surface)
{
  float x_orig = x;
  float y_orig = y;
  x -= x_idx * cell_width;
  y -= y_idx * cell_width;

  int current_surface = 0;

  RT_FLOAT min_dist = 1e9;
  RT_FLOAT surface_normal_x = 0.0;
  RT_FLOAT surface_normal_y = 0.0;

  // Test out all 4 surfaces
  RT_FLOAT x_pos_dist = (cell_width - x) / x_dir;
  RT_FLOAT y_pos_dist = (cell_width - y) / y_dir;
  RT_FLOAT x_neg_dist = -x / x_dir;
  RT_FLOAT y_neg_dist = -y / y_dir;

  // Determine closest one
  if( x_pos_dist < min_dist && x_pos_dist > 0 && last_surface != -1 )
  {
    min_dist = x_pos_dist;
    surface_normal_x = 1.0;
    surface_normal_y = 0.0;
    current_surface = 1;
  }
  if( y_pos_dist < min_dist && y_pos_dist > 0 && last_surface != -2)
  {
    min_dist = y_pos_dist;
    surface_normal_x = 0.0;
    surface_normal_y = 1.0;
    current_surface = 2;
  }
  if( x_neg_dist < min_dist && x_neg_dist > 0 && last_surface != 1)
  {
    min_dist = x_neg_dist;
    surface_normal_x = -1.0;
    surface_normal_y =  0.0;
    current_surface = -1;
  }
  if( y_neg_dist < min_dist && y_neg_dist > 0 && last_surface != 2)
  {
    min_dist = y_neg_dist;
    surface_normal_x =  0.0;
    surface_normal_y = -1.0;
    current_surface = -2;
  }

  assert(min_dist > 0 );
  if( min_dist > 1e2)
  {
    printf("[x+, y+, x-, y-] = [%e, %e, %e, %e] last surface = %d  x_relative = %e   y_relative = %e   x_orig = %e   y_orig = %e  x_idx = %d  y_idx = %d   cell_width = %e \n", x_pos_dist, y_pos_dist, x_neg_dist, y_neg_dist, last_surface, x, y, x_orig, y_orig, x_idx, y_idx, cell_width);
    printf("min dist = %f\n", min_dist);
  }
  assert(min_dist < 1e2);
  TraceResult trace;
  trace.distance_to_surface = min_dist;
  trace.surface_normal_x = surface_normal_x;
  trace.surface_normal_y = surface_normal_y;
  trace.last_surface = current_surface;

  return trace;
}
