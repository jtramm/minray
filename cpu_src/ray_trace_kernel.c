#include "minray.h"

void ray_trace_kernel(Parameters P, SimulationData SD, RayData rayData, uint64_t ray_id)
{
  double distance_travelled = 0.0;
  int intersection_id = 0;

  double x =     rayData.location_x[ ray_id];
  double y =     rayData.location_y[ ray_id];
  double x_dir = rayData.direction_x[ray_id];
  double y_dir = rayData.direction_y[ray_id];
  int cell_id =  rayData.cell_id[    ray_id];
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
    TraceResult trace = cartesian_ray_trace(x, y, P.cell_width, x_idx, y_idx, x_dir, y_dir);
  
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
    double x_across_surface = x + trace.surface_normal_x * BUMP;
    double y_across_surface = y + trace.surface_normal_y * BUMP;
    
    // Look up the "neighbor" cell id of the test point.
    // This function also gives us some info on if we hit a boundary, and what type it was.
    //CellLookup lookup = find_cell_id(P, x_across_surface, y_across_surface);
    CellLookup lookup = find_cell_id_using_neighbor_list(P, &SD.readWriteData.cellData.neighborList[cell_id], x_across_surface, y_across_surface);
    
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
    x += trace.surface_normal_x * BUMP;
    y += trace.surface_normal_y * BUMP;

    // Add this intersection's distance to the total for the ray
    distance_travelled += trace.distance_to_surface;

    // Some sanity checks (can be disabled if desired)
    assert(cell_id >= 0 && cell_id < P.n_cells);
    assert(x > 0.0 && y > 0.0 && x < P.length_per_dimension && y < P.length_per_dimension);
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
    
  // Bank number of intersections that this ray had this iteration
  SD.readWriteData.intersectionData.n_intersections[ray_id] = intersection_id;
}

int is_point_inside_CSG_cell(Parameters P, double x, double y, int cell_id)
{
  int cartesian_cell_idx_x = floor(x * P.inverse_cell_width);
  int cartesian_cell_idx_y = floor(y * P.inverse_cell_width);
  int true_cell_id = cartesian_cell_idx_y * P.n_cells_per_dimension + cartesian_cell_idx_x; 
  if( cell_id == true_cell_id )
    return 1;
  else
    return 0;
}

CellLookup find_cell_id_general_fast(Parameters P, double x, double y)
{
  int cartesian_cell_idx_x = floor(x * P.inverse_cell_width);
  int cartesian_cell_idx_y = floor(y * P.inverse_cell_width);

  int boundary_x = floor(x * P.inverse_length_per_dimension) + 1;
  int boundary_y = floor(y * P.inverse_length_per_dimension) + 1;

  int boundary_condition = P.boundary_conditions[boundary_x][boundary_y];

  int cell_id = cartesian_cell_idx_y * P.n_cells_per_dimension + cartesian_cell_idx_x; 

  CellLookup lookup;
  lookup.cell_id = cell_id;
  lookup.cartesian_cell_idx_x = cartesian_cell_idx_x;
  lookup.cartesian_cell_idx_y = cartesian_cell_idx_y;
  lookup.boundary_condition = boundary_condition;
  return lookup;
}

CellLookup find_cell_id_general(Parameters P, double x, double y)
{
  int boundary_x = floor(x * P.inverse_length_per_dimension) + 1;
  int boundary_y = floor(y * P.inverse_length_per_dimension) + 1;
  int boundary_condition = P.boundary_conditions[boundary_x][boundary_y];

  int cell_id = -1;
  for( int i = 0; i < P.n_cells; i++ )
  {
    if( is_point_inside_CSG_cell(P, x, y, i) )
    {
      cell_id = i;
      break;
    }
  }

  CellLookup lookup;
  lookup.cell_id = cell_id;
  lookup.cartesian_cell_idx_x = cell_id % P.n_cells_per_dimension;
  lookup.cartesian_cell_idx_y = cell_id / P.n_cells_per_dimension;
  lookup.boundary_condition = boundary_condition;
  return lookup;
}


CellLookup find_cell_id_using_neighbor_list(Parameters P, NeighborList * neighborList, double x, double y)
{
  // Determine boundary information
  int boundary_x = floor(x * P.inverse_length_per_dimension) + 1;
  int boundary_y = floor(y * P.inverse_length_per_dimension) + 1;
  int boundary_condition = P.boundary_conditions[boundary_x][boundary_y];
  if( boundary_condition != NONE )
  {
    CellLookup lookup;
    lookup.cell_id = -1;
    lookup.boundary_condition = boundary_condition;
    return lookup;
  }

  NeighborListIterator iterator;
  nl_init_iterator(neighborList, &iterator);

  int cell_id = -1;

  // Test (x,y) location against all CSG cells in the neighbor list
  while(!iterator.is_finished)
  {
    int neighbor_id = nl_read_next(neighborList, &iterator);
    assert(neighbor_id >= 0 && neighbor_id < P.n_cells);

    if( is_point_inside_CSG_cell(P, x, y, neighbor_id) )
    {
      cell_id = neighbor_id;
      break;
    }
  }

  // If (x,y) location is not found in neighbor list, perform general search and add result to neighbor list
  if( cell_id == -1 )
  {
    CellLookup lookup = find_cell_id_general(P, x, y);
    nl_push_back(neighborList, lookup.cell_id);
    return lookup;
  }

  CellLookup lookup;
  lookup.cell_id = cell_id;
  lookup.cartesian_cell_idx_x = cell_id % P.n_cells_per_dimension;
  lookup.cartesian_cell_idx_y = cell_id / P.n_cells_per_dimension;
  lookup.boundary_condition = boundary_condition;
  return lookup;
}

TraceResult cartesian_ray_trace(double x, double y, double cell_width, int x_idx, int y_idx, double x_dir, double y_dir)
{
  x -= x_idx * cell_width;
  y -= y_idx * cell_width;

  double min_dist = 1e9;
  double surface_normal_x = 0.0;
  double surface_normal_y = 0.0;

  // Test out all 4 surfaces
  double x_pos_dist = (cell_width - x) / x_dir;
  double y_pos_dist = (cell_width - y) / y_dir;
  double x_neg_dist = -x / x_dir;
  double y_neg_dist = -y / y_dir;

  // Determine closest one
  if( x_pos_dist < min_dist && x_pos_dist > 0 )
  {
    min_dist = x_pos_dist;
    surface_normal_x = 1.0;
    surface_normal_y = 0.0;
  }
  if( y_pos_dist < min_dist && y_pos_dist > 0 )
  {
    min_dist = y_pos_dist;
    surface_normal_x = 0.0;
    surface_normal_y = 1.0;
  }
  if( x_neg_dist < min_dist && x_neg_dist > 0)
  {
    min_dist = x_neg_dist;
    surface_normal_x = -1.0;
    surface_normal_y =  0.0;
  }
  if( y_neg_dist < min_dist && y_neg_dist > 0)
  {
    min_dist = y_neg_dist;
    surface_normal_x =  0.0;
    surface_normal_y = -1.0;
  }

  TraceResult trace;
  trace.distance_to_surface = min_dist;
  trace.surface_normal_x = surface_normal_x;
  trace.surface_normal_y = surface_normal_y;

  return trace;
}
