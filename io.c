#include "minray.h"

Parameters read_CLI(int argc, char * argv[])
{
  Parameters P;
  // User Inputs
  P.n_cells_per_dimension = 102;
  P.length_per_dimension = 64.26;
  //P.height = 100.0;
  P.n_rays = 7500;
  P.distance_per_ray = 10.0;
  P.n_inactive_iterations = 100;
  P.n_active_iterations = 100;
  P.seed = 1337;
  P.n_materials = 8;
  P.n_energy_groups = 7;
  P.max_intersections_per_ray = 1000;
  P.boundary_conditions[0] = NONE;
  P.boundary_conditions[1] = REFLECTIVE;
  P.boundary_conditions[2] = VACUUM;
  P.boundary_conditions[3] = REFLECTIVE;
  P.boundary_conditions[4] = VACUUM;

  // Derived Values
  P.cell_width = P.length_per_dimension / P.n_cells_per_dimension;
  P.inverse_cell_width = 1.0 / P.cell_width;
  P.n_cells = P.n_cells_per_dimension * P.n_cells_per_dimension;
  P.cell_expected_track_length = (P.distance_per_ray * P.n_rays) / P.n_cells;
  P.inverse_total_track_length = 1.0 / (P.distance_per_ray * P.n_rays);
  P.inverse_length_per_dimension = 1.0 / P.length_per_dimension;
  P.n_iterations = P.n_inactive_iterations + P.n_active_iterations;

  return P;
}

// KEY
// 0 - UO2
// 1 - MOX 4.3
// 2 - MOX 7.0
// 3 - MOX 8.7
// 4 - Fission Chamber
// 5 - Guide Tube
// 6 - Moderator 
// 7 - Control Rod
ReadOnlyData load_2D_C5G7_XS(Parameters P)
{
  size_t sz = P.n_materials * P.n_energy_groups;

  // Allocate Data on Host
  float * nu_Sigma_f  = (float *) calloc(sz, sizeof(float));
  float * Sigma_f     = (float *) calloc(sz, sizeof(float));
  float * Sigma_t     = (float *) calloc(sz, sizeof(float));
  float * Chi         = (float *) calloc(sz, sizeof(float));

  sz = sz * P.n_energy_groups;
  float * Sigma_s =     (float *) calloc(sz, sizeof(float));

  // UO2 fuel-clad mixture
  FILE * UO2_scatter   = fopen("data/C5G7_2D/UO2_scatter.txt",   "r");
  FILE * UO2_transport = fopen("data/C5G7_2D/UO2_transport.txt", "r");
  FILE * UO2_nufission = fopen("data/C5G7_2D/UO2_nufission.txt", "r");
  FILE * UO2_fission   = fopen("data/C5G7_2D/UO2_fission.txt", "r");
  FILE * UO2_chi       = fopen("data/C5G7_2D/UO2_chi.txt", "r");

  // MOX 4.3% fuel-clad mixture
  FILE * MOX_43_scatter   = fopen("data/C5G7_2D/MOX_43_scatter.txt",   "r");
  FILE * MOX_43_transport = fopen("data/C5G7_2D/MOX_43_transport.txt", "r");
  FILE * MOX_43_nufission = fopen("data/C5G7_2D/MOX_43_nufission.txt", "r");
  FILE * MOX_43_fission   = fopen("data/C5G7_2D/MOX_43_fission.txt", "r");
  FILE * MOX_43_chi       = fopen("data/C5G7_2D/MOX_43_chi.txt", "r");

  // MOX 7.0% fuel-clad mixture
  FILE * MOX_70_scatter   = fopen("data/C5G7_2D/MOX_70_scatter.txt",   "r");
  FILE * MOX_70_transport = fopen("data/C5G7_2D/MOX_70_transport.txt", "r");
  FILE * MOX_70_nufission = fopen("data/C5G7_2D/MOX_70_nufission.txt", "r");
  FILE * MOX_70_fission   = fopen("data/C5G7_2D/MOX_70_fission.txt", "r");
  FILE * MOX_70_chi       = fopen("data/C5G7_2D/MOX_70_chi.txt", "r");

  // MOX 8.7% fuel-clad mixture
  FILE * MOX_87_scatter   = fopen("data/C5G7_2D/MOX_87_scatter.txt",   "r");
  FILE * MOX_87_transport = fopen("data/C5G7_2D/MOX_87_transport.txt", "r");
  FILE * MOX_87_nufission = fopen("data/C5G7_2D/MOX_87_nufission.txt", "r");
  FILE * MOX_87_fission   = fopen("data/C5G7_2D/MOX_87_fission.txt", "r");
  FILE * MOX_87_chi       = fopen("data/C5G7_2D/MOX_87_chi.txt", "r");

  // Fission Chamber
  FILE * FC_scatter   = fopen("data/C5G7_2D/FC_scatter.txt",   "r");
  FILE * FC_transport = fopen("data/C5G7_2D/FC_transport.txt", "r");
  FILE * FC_nufission = fopen("data/C5G7_2D/FC_nufission.txt", "r");
  FILE * FC_fission   = fopen("data/C5G7_2D/FC_fission.txt", "r");
  FILE * FC_chi       = fopen("data/C5G7_2D/FC_chi.txt", "r");

  // Guide Tube
  FILE * GT_scatter   = fopen("data/C5G7_2D/GT_scatter.txt",   "r");
  FILE * GT_transport = fopen("data/C5G7_2D/GT_transport.txt", "r");

  // Moderator
  FILE * Mod_scatter   = fopen("data/C5G7_2D/Mod_scatter.txt",   "r");
  FILE * Mod_transport = fopen("data/C5G7_2D/Mod_transport.txt", "r");

  // Control Rod
  FILE * CR_scatter   = fopen("data/C5G7_2D/CR_scatter.txt",   "r");
  FILE * CR_transport = fopen("data/C5G7_2D/CR_transport.txt", "r");

  int ret = 0;

  for( int i = 0; i < P.n_energy_groups; i++ )
  {
    // Transport XS
    ret = fscanf(UO2_transport,    "%f", Sigma_t + 0 * P.n_energy_groups + i);
    ret = fscanf(MOX_43_transport, "%f", Sigma_t + 1 * P.n_energy_groups + i);
    ret = fscanf(MOX_70_transport, "%f", Sigma_t + 2 * P.n_energy_groups + i);
    ret = fscanf(MOX_87_transport, "%f", Sigma_t + 3 * P.n_energy_groups + i);
    ret = fscanf(FC_transport,     "%f", Sigma_t + 4 * P.n_energy_groups + i);
    ret = fscanf(GT_transport,     "%f", Sigma_t + 5 * P.n_energy_groups + i);
    ret = fscanf(Mod_transport,    "%f", Sigma_t + 6 * P.n_energy_groups + i);
    ret = fscanf(CR_transport,     "%f", Sigma_t + 7 * P.n_energy_groups + i);

    // Nu-Fission XS
    ret = fscanf(UO2_nufission,    "%f", nu_Sigma_f + 0 * P.n_energy_groups + i);
    ret = fscanf(MOX_43_nufission, "%f", nu_Sigma_f + 1 * P.n_energy_groups + i);
    ret = fscanf(MOX_70_nufission, "%f", nu_Sigma_f + 2 * P.n_energy_groups + i);
    ret = fscanf(MOX_87_nufission, "%f", nu_Sigma_f + 3 * P.n_energy_groups + i);
    ret = fscanf(FC_nufission,     "%f", nu_Sigma_f + 4 * P.n_energy_groups + i);

    // Fission XS
    ret = fscanf(UO2_fission,    "%f", Sigma_f + 0 * P.n_energy_groups + i);
    ret = fscanf(MOX_43_fission, "%f", Sigma_f + 1 * P.n_energy_groups + i);
    ret = fscanf(MOX_70_fission, "%f", Sigma_f + 2 * P.n_energy_groups + i);
    ret = fscanf(MOX_87_fission, "%f", Sigma_f + 3 * P.n_energy_groups + i);
    ret = fscanf(FC_fission,     "%f", Sigma_f + 4 * P.n_energy_groups + i);

    // Chi
    ret = fscanf(UO2_chi,    "%f", Chi + 0 * P.n_energy_groups + i);
    ret = fscanf(MOX_43_chi, "%f", Chi + 1 * P.n_energy_groups + i);
    ret = fscanf(MOX_70_chi, "%f", Chi + 2 * P.n_energy_groups + i);
    ret = fscanf(MOX_87_chi, "%f", Chi + 3 * P.n_energy_groups + i);
    ret = fscanf(FC_chi,     "%f", Chi + 4 * P.n_energy_groups + i);

    for( int j = 0; j < P.n_energy_groups; j++ )
    {
      // Scatter XS
      ret = fscanf(UO2_scatter,    "%f", Sigma_s + 0 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(MOX_43_scatter, "%f", Sigma_s + 1 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(MOX_70_scatter, "%f", Sigma_s + 2 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(MOX_87_scatter, "%f", Sigma_s + 3 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(FC_scatter,     "%f", Sigma_s + 4 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(GT_scatter,     "%f", Sigma_s + 5 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(Mod_scatter,    "%f", Sigma_s + 6 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
      ret = fscanf(CR_scatter,     "%f", Sigma_s + 7 * P.n_energy_groups * P.n_energy_groups + i * P.n_energy_groups + j);
    }
  }

  // Close all Files
  fclose(UO2_scatter      );
  fclose(UO2_transport    );
  fclose(UO2_nufission    );
  fclose(UO2_fission    );
  fclose(UO2_chi          );
  fclose(MOX_43_scatter   );
  fclose(MOX_43_transport );
  fclose(MOX_43_nufission );
  fclose(MOX_43_fission );
  fclose(MOX_43_chi       );
  fclose(MOX_70_scatter   );
  fclose(MOX_70_transport );
  fclose(MOX_70_nufission );
  fclose(MOX_70_fission );
  fclose(MOX_70_chi       );
  fclose(MOX_87_scatter   );
  fclose(MOX_87_transport );
  fclose(MOX_87_nufission );
  fclose(MOX_87_fission );
  fclose(MOX_87_chi       );
  fclose(FC_scatter       );
  fclose(FC_transport     );
  fclose(FC_nufission     );
  fclose(FC_fission     );
  fclose(FC_chi           );
  fclose(GT_scatter       );
  fclose(GT_transport     );
  fclose(Mod_scatter      );
  fclose(Mod_transport    );
  fclose(CR_scatter      );
  fclose(CR_transport    );

  char fname[512];
  sprintf(fname, "data/C5G7_2D/material_ids_%d.txt", P.n_cells_per_dimension);
  FILE * material_file = fopen(fname, "r");
  sz = P.n_cells * sizeof(int);
  int * material_id = (int *) malloc(sz);
  for( int c = 0; c < P.n_cells; c++ )
  {
    ret = fscanf(material_file, "%d", material_id + c);
  }

  fclose(material_file);

  ReadOnlyData ROD;
  ROD.Sigma_f = Sigma_f;
  ROD.Sigma_t = Sigma_t;
  ROD.Sigma_s = Sigma_s;
  ROD.nu_Sigma_f = nu_Sigma_f;
  ROD.Chi = Chi;
  ROD.material_id = material_id;

  if( ret == 0 )
  {
    printf("something went wrong with XS read in...\n");
    exit(1);
  }

  return ROD;
}

float eswap_float( float f )
{    
  char * ptr = (char *) &f;
  char orig[4];
  orig[0] = ptr[0];
  orig[1] = ptr[1];
  orig[2] = ptr[2];
  orig[3] = ptr[3];
  char new[4];
  new[0] = orig[3];
  new[1] = orig[2];
  new[2] = orig[1];
  new[3] = orig[0];
  float f2 = *((float *)new);
  return f2;
}

int eswap_int( int f )
{    
  char * ptr = (char *) &f;
  char orig[4];
  orig[0] = ptr[0];
  orig[1] = ptr[1];
  orig[2] = ptr[2];
  orig[3] = ptr[3];
  char new[4];
  new[0] = orig[3];
  new[1] = orig[2];
  new[2] = orig[1];
  new[3] = orig[0];
  int f2 = *((int *)new);
  return f2;
}

void plot_3D_vtk(Parameters P, float * scalar_flux_accumulator, int * material_id)
{
  int N = P.n_cells_per_dimension;
  int z_N = 1;

  char fname[512];
  sprintf(fname, "minray_plots_%d.vtk", N);
  FILE * fp = fopen(fname, "w");

  int plot_thermal_flux = 1;
  int plot_fast_flux = 1;
  int plot_materials = 1;

  /*
  // For Corner Zoom
  G.box.min.x = 0.0 + 17*1.26 / 2.0 - 2*1.26;
  G.box.max.x = G.box.min.x + 4 * 1.26;
  G.box.min.y = 0.0 - 17*1.26 / 2.0 - 2*1.26;
  G.box.max.y = G.box.min.y + 4 * 1.26;
  */

  double x_delta = P.cell_width;
  double y_delta = P.cell_width;
  double z_delta = P.cell_width;

  printf("plotting 3D Data X x Y x Z = %d x %d x %d\n", N, N, z_N);

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp, "Dataset File\n");
  //fprintf(fp, "ASCII\n");
  fprintf(fp, "BINARY\n");
  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %d %d %d\n", N, N, z_N);
  fprintf(fp, "ORIGIN 0 0 0\n");
  fprintf(fp, "SPACING %lf %lf %lf\n", x_delta, y_delta, z_delta);
  fprintf(fp, "POINT_DATA %d\n", N*N*z_N);

  if( plot_thermal_flux )
  {
    fprintf(fp, "SCALARS thermal_flux float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    int cell_id = 0;
    for( int y = 0; y < P.n_cells_per_dimension; y++)
    {
      for( int x = 0; x < P.n_cells_per_dimension; x++)
      {
        float thermal_flux = scalar_flux_accumulator[cell_id * P.n_energy_groups] / P.n_active_iterations;
        thermal_flux = eswap_float(thermal_flux);
        fwrite(&thermal_flux, sizeof(float), 1, fp);
        cell_id++;
      }
    }
  }
  
  if( plot_fast_flux )
  {
    fprintf(fp, "SCALARS fast_flux float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    int cell_id = 0;
    for( int y = 0; y < P.n_cells_per_dimension; y++)
    {
      for( int x = 0; x < P.n_cells_per_dimension; x++)
      {
        float fast_flux = scalar_flux_accumulator[cell_id * P.n_energy_groups + P.n_energy_groups - 1] / P.n_active_iterations;
        fast_flux = eswap_float(fast_flux);
        fwrite(&fast_flux, sizeof(float), 1, fp);
        cell_id++;
      }
    }
  }

  if( plot_materials )
  {
    fprintf(fp, "SCALARS material_type int\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    int cell_id = 0;
    for( int y = 0; y < P.n_cells_per_dimension; y++)
    {
      for( int x = 0; x < P.n_cells_per_dimension; x++)
      {
        int material = material_id[cell_id++];
        material = eswap_int(material);
        fwrite(&material, sizeof(int), 1, fp);
      }
    }
  }

  fclose(fp);
  printf("Finished plotting!\n");
}

