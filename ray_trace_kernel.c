#include "minray.h"

//#define BUMP 1.0e-12
//#define BUMP 1.0e-11
#define BUMP 1.0e-10
//#define BUMP 1.0e-9

int find_cell_id(double x, double y, double inverse_cell_width, int n_cells_per_dimension, double inverse_length_per_dimension, int * boundary_surface, int * x_idx, int * y_idx);
double cartesian_ray_trace(double x, double y, double cell_width, int x_idx, int y_idx, double x_dir, double y_dir, int * intersected_surface_direction);

void ray_trace_kernel(Parameters P, SimulationData SD, RayData rayData, uint64_t ray_id)
{
  double distance_travelled = 0.0;
  int intersection_id = 0;
  int boundary_surface = 0;

  double x =     rayData.location_x[ray_id];
  double y =     rayData.location_y[ray_id];
  double x_dir = rayData.direction_x[ray_id];
  double y_dir = rayData.direction_y[ray_id];
  //double z_dir = rayData.direction_z[ray_id];
  int cell_id =  rayData.cell_id[ray_id];
  int x_idx = cell_id % P.n_cells_per_dimension;
  int y_idx = cell_id / P.n_cells_per_dimension;

  int just_hit_vacuum = 0;
  int is_terminal = 0;

  // We run this loop until either:
  // 1) The maximum number of intersections has been reached
  // 2) The ray has reached its set distance
  for( intersection_id = 0; (intersection_id < P.max_intersections_per_ray) && (distance_travelled < P.distance_per_ray); intersection_id++ )
  {
    //printf("beginning intersection %d\n", intersection_id);
    //print_ray(x, y, x_dir, y_dir, cell_id);

    // Find distance to next surface
    int intersected_surface_direction = 0;
    double distance_to_surface = cartesian_ray_trace(x, y, P.cell_width, x_idx, y_idx, x_dir, y_dir, &intersected_surface_direction);
    //printf("distance to surface = %.3lf\n", distance_to_surface);
  
    // Check to see if we are terminal
    if(distance_travelled + distance_to_surface >= P.distance_per_ray)
    {
      distance_to_surface = (P.distance_per_ray - distance_travelled) + BUMP;
      is_terminal = 1;
    }
  
    // Record Intersection information
    uint64_t global_intersection_id = ray_id * P.max_intersections_per_ray + intersection_id;

    // Record distance to next surface
    SD.readWriteData.intersectionData.distances[global_intersection_id] = distance_to_surface;

    // Record the ID of the cell we are in
    SD.readWriteData.intersectionData.cell_ids[global_intersection_id] = cell_id;
    
    // Record if we hit vacuum
    SD.readWriteData.intersectionData.did_vacuum_reflects[global_intersection_id] = just_hit_vacuum;
    just_hit_vacuum = 0;
    
    // Record if we hit vacuum
    // TODO NOTE: If running in parallel, this may need to be an atomic operation.
    SD.readWriteData.cellData.hit_count[cell_id] = 1;

    // Move ray forward to intersection surface
    x += x_dir * distance_to_surface;
    y += y_dir * distance_to_surface;

    /*
    double x_across_surface = x + x_dir * BUMP;
    double y_across_surface = y + y_dir * BUMP;
    */
    double x_across_surface = x;
    double y_across_surface = y;
    if( intersected_surface_direction == X_NEG )
      x_across_surface -= BUMP;
    else if( intersected_surface_direction == X_POS )
      x_across_surface += BUMP;
    else if( intersected_surface_direction == Y_NEG )
      y_across_surface -= BUMP;
    else if( intersected_surface_direction == Y_POS )
      y_across_surface += BUMP;
    int x_idx_across_surface, y_idx_across_surface;
    
    // Look up cell_ID of neighbor we are travelling into
    int neighbor_id = find_cell_id(x_across_surface, y_across_surface, P.inverse_cell_width, P.n_cells_per_dimension, P.inverse_length_per_dimension, &boundary_surface, &x_idx_across_surface, &y_idx_across_surface);
    //printf("neighbor_id = %d\n", neighbor_id);
    if( neighbor_id == cell_id && !is_terminal)
    {
      printf("neighbor_id = cell_id = %d (x_idx_ac = %d y_idx_ac = %d)\n", neighbor_id, x_idx_across_surface, y_idx_across_surface);
      printf("neighbor_id = cell_id = %d (x_idx    = %d y_idx    = %d)\n", neighbor_id, x_idx, y_idx);
      printf("x, y = [%le, %le]  neib_x, neib_y = [%le, %le])\n", x, y, x_across_surface, y_across_surface);
      printf("Location within cell = [%lf, %lf]\n", x - x_idx_across_surface* P.cell_width, y - y_idx_across_surface * P.cell_width);
      printf("distance travelled to surface = %lf\n", distance_to_surface);
      printf("intersection_id = %d\n", intersection_id);
      print_ray(x, y, x_dir, y_dir, cell_id);
      assert(neighbor_id != cell_id);
    }

    
    //printf("boundary surface = %d\n", boundary_surface);
    // Reflect
    if( boundary_surface )
    {
      //printf("just hit a boundary\n");
      if( boundary_surface % 2 == 1 )
        x_dir *= -1.0;
      else
        y_dir *= -1.0;
    }

    // Lookup boundary condition information based on boundary surface hit
    int boundary_condition = P.boundary_conditions[boundary_surface];

    // Note boundary information for next intersection.
    if( boundary_condition == VACUUM )
    {
      just_hit_vacuum = 1;
      //printf("Just hit vacuum boundary\n");
    }

    // If we didn't hit a boundary, move into the next cell.
    if( !boundary_condition )
    {
      cell_id = neighbor_id;
      x_idx = x_idx_across_surface;
      y_idx = y_idx_across_surface;
    }
    
    assert(cell_id >= 0 && cell_id < P.n_cells);

    // Move ray off the surface
    /*
    x += x_dir * BUMP;
    y += y_dir * BUMP;
    */
    
    if( boundary_condition == NONE )
    {
      if( intersected_surface_direction == X_NEG )
        x -= BUMP;
      else if( intersected_surface_direction == X_POS )
        x += BUMP;
      else if( intersected_surface_direction == Y_NEG )
        y -= BUMP;
      else if( intersected_surface_direction == Y_POS )
        y += BUMP;
    }
    else
    {
      if( intersected_surface_direction == X_NEG )
        x += BUMP;
      else if( intersected_surface_direction == X_POS )
        x -= BUMP;
      else if( intersected_surface_direction == Y_NEG )
        y += BUMP;
      else if( intersected_surface_direction == Y_POS )
        y -= BUMP;
    }
    
    assert(x > 0.0 && y > 0.0 && x < P.length_per_dimension && y < P.length_per_dimension);
    //printf("Location within cell = [%lf, %lf]\n", x - x_idx_across_surface* P.cell_width, y - y_idx_across_surface * P.cell_width);


    distance_travelled += distance_to_surface;
    //printf("distance_travelled = %lf\n", distance_travelled);
    //print_ray(x, y, x_dir, y_dir, cell_id);
  }
  //assert(intersection_id < P.max_intersections_per_ray);
  if(intersection_id >= P.max_intersections_per_ray)
  {
    printf("WARNING: Increase max number of intersections per ray\n");
    print_ray(x, y, x_dir, y_dir, cell_id);
  }
  
  
  rayData.location_x[ray_id] = x;
  rayData.location_y[ray_id] = y;
  rayData.direction_x[ray_id] = x_dir;
  rayData.direction_y[ray_id] = y_dir;
  //rayData.direction_z[ray_id] = z_dir;
  rayData.cell_id[ray_id] = cell_id;
    
  SD.readWriteData.intersectionData.n_intersections[ray_id] = intersection_id;
}

// Boundary surfaces:
// inside_domain = 0
// negative_x = 1;
// positive_x = 3;
// negative_y = 2;
// positive_y = 4;
int find_cell_id(double x, double y, double inverse_cell_width, int n_cells_per_dimension, double inverse_length_per_dimension, int * boundary_surface, int * x_idx, int * y_idx)
{
  *x_idx = floor(x * inverse_cell_width);
  *y_idx = floor(y * inverse_cell_width);
  //printf("x_idx = %d, y_idx = %d\n", *x_idx, *y_idx);

  int boundary_x = floor(x * inverse_length_per_dimension);
  int boundary_y = floor(y * inverse_length_per_dimension);
  //printf("[%lf, %lf]  boundary_x = %d boundary_y = %d inverse_length_per = %lf\n", x, y, boundary_x, boundary_y, inverse_length_per_dimension);

  if( boundary_x )
    *boundary_surface = boundary_x+2;
  else if( boundary_y )
    *boundary_surface = boundary_y+3;
  else
    *boundary_surface = 0;

  int cell_id = *y_idx * n_cells_per_dimension + *x_idx; 
  return cell_id;
}

double cartesian_ray_trace(double x, double y, double cell_width, int x_idx, int y_idx, double x_dir, double y_dir, int * intersected_surface_direction)
{
  x -= x_idx * cell_width;
  y -= y_idx * cell_width;

  double min_dist = 1e9;
  //printf("Shifted coord [%lf, %lf]\n", x, y);
  //printf("X,Y Indices:  [%d, %d]\n", x_idx, y_idx);

  // Test out all 4 surfaces
  double x_pos_dist = (cell_width - x) / x_dir;
  double y_pos_dist = (cell_width - y) / y_dir;
  double x_neg_dist = -x / x_dir;
  double y_neg_dist = -y / y_dir;

  // Determine closest one
  if( x_pos_dist < min_dist && x_pos_dist > 0 )
  {
    min_dist = x_pos_dist;
    *intersected_surface_direction = X_POS;
  }
  if( y_pos_dist < min_dist && y_pos_dist > 0 )
  {
    min_dist = y_pos_dist;
    *intersected_surface_direction = Y_POS;
  }
  if( x_neg_dist < min_dist && x_neg_dist > 0)
  {
    min_dist = x_neg_dist;
    *intersected_surface_direction = X_NEG;
  }
  if( y_neg_dist < min_dist && y_neg_dist > 0)
  {
    min_dist = y_neg_dist;
    *intersected_surface_direction = Y_NEG;
  }

  return min_dist;
}

