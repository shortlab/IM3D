/******************************************************************************
  Module Name : target.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Contains the target-related functions, etc.

  Others :
      Error codes used in this modul: between 3000 and 3999.
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "target.h"

/*=============================================================================
  Function Name : init_target_structure
  Description   : Read target structure (size etc.) from config file and reads
                  target concentration array etc.
                  First, we need to read the config file about how large the
                  target is, how many cells and so on... .
                  Then we can initialze the target arrays, and then finally read
                  in the composition matrix.

  Inputs  : char* file_name
  Outputs : no.

  Notes :
      TODO: free arrays
      Function call:
        fileio - int read_init_file (int (*read_data_block) (char* block_name),
                            int (*read_data) (char* par_name, char* par_value),
                            char* file_name);
                 int read_float_file_into_array (char* file_name, float* target_array,
                            int count, int file_type)
        util - void get_float_one_bit_smaller (float* flt_input, float* flt_output);
        csg  - int read_csg_shape ();
        fetm - int read_fetm_shape ().
=============================================================================*/
int init_target_structure (char *file_name) {
    int i, j, result;
#ifndef MPI_PRALLEL
    /* MPI=============================================== */
    unsigned long int lui_temp = 0;
    /* MPI=============================================== */
#endif

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=0) printf ("Initializing target.\n");
    /* MPI=============================================== */
#else
    if (print_level >= 0) printf ("Initializing target.\n");
#endif

    /* set some default values */
    use_density_mult = 0;

    /* read structure definition file */  /* (1) - (8) - (1) */
    result = read_init_file (read_target_structure_data_block, read_target_structure_data, file_name);
    if (result != 0) return result;
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=0) printf ("Target structure definition file: %s\n", file_name);
    /* MPI=============================================== */
#else
    if (print_level >= 0) printf ("Target structure definition file: %s\n", file_name);
#endif

    /* Right hand coordinate system: x-axis vertical to target surface (point into substrate),
       y,z-axises parallel to target surface. Uniform cell is only for the special target, for
       the substrate regions outside of the target region, nonuniform cells will be used or do
       not be counted. Or, this cell group is indepedent of the geomentry cell group, so the
       range of this cell group can be changed arbitrarily. */
    layer_count_yz = cell_count_y   * cell_count_z;
    cell_count     = layer_count_yz * cell_count_x;
    target_size_x  = cell_count_x   * cell_size_x;
    target_size_y  = cell_count_y   * cell_size_y;
    target_size_z  = cell_count_z   * cell_size_z;

    /* for nodes in msh module */
    node_x = cell_count_x + 1;
    node_y = cell_count_y + 1;
    node_z = cell_count_z + 1;
    node_yz = node_y * node_z;
    node_count = node_x * node_yz;

    /* we sometimes need a float that is one bit smaller than the targetsize itself */
    get_float_one_bit_smaller (&target_size_x, &target_max_x);
    get_float_one_bit_smaller (&target_size_y, &target_max_y);
    get_float_one_bit_smaller (&target_size_z, &target_max_z);

    cell_volume = cell_size_x * cell_size_y * cell_size_z;
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=0) {
        printf ("\nTarget size is: \n");
        printf ("x: %i cells, %g nm per cell, %g nm in total.\n",
                cell_count_x, cell_size_x, target_size_x);
        printf ("y: %i cells, %g nm per cell, %g nm in total.\n",
                cell_count_y, cell_size_y, target_size_y);
        printf ("z: %i cells, %g nm per cell, %g nm in total.\n",
                cell_count_z, cell_size_z, target_size_z);
        printf ("Total: %i cells in %g nm^3.\n", cell_count,
                target_size_x * target_size_y * target_size_z);
    }
    /* MPI=============================================== */
#else
    if (print_level >= 0) {
        printf ("\nTarget size is: \n");
        printf ("x: %i cells, %g nm per cell, %g nm in total.\n",
                cell_count_x, cell_size_x, target_size_x);
        printf ("y: %i cells, %g nm per cell, %g nm in total.\n",
                cell_count_y, cell_size_y, target_size_y);
        printf ("z: %i cells, %g nm per cell, %g nm in total.\n",
                cell_count_z, cell_size_z, target_size_z);
        printf ("Total: %i cells in %g nm^3.\n", cell_count,
                target_size_x * target_size_y * target_size_z);
    }
#endif

    /* Now that we know the size of the target we can initialze some arrays etc. */
    if (mem_usage_only == 0) {
        /* Use calloc instead of malloc for initializing to zeros */
        target_composition  = (int*) calloc (cell_count, sizeof (int));
        if (target_composition  == NULL) return -3001;  /* cannot allocate memory */
        target_implanted_ions = (int*) calloc (cell_count, sizeof (int));
        target_replacing_ions = (int*) calloc (cell_count, sizeof (int));
        target_total_displacements = (int*) calloc (cell_count, sizeof (int));
        target_total_interstitials = (int*) calloc (cell_count, sizeof (int));
        target_total_vacancies     = (int*) calloc (cell_count, sizeof (int));
        target_total_replacements  = (int*) calloc (cell_count, sizeof (int));

        target_depth_implanted_ions = (int*) calloc (cell_count_z, sizeof (int));;
        target_depth_replacing_ions = (int*) calloc (cell_count_z, sizeof (int));
        target_depth_total_displacements = (int*) calloc (cell_count_z, sizeof (int));
        target_depth_total_interstitials = (int*) calloc (cell_count_z, sizeof (int));
        target_depth_total_vacancies     = (int*) calloc (cell_count_z, sizeof (int));
        target_depth_total_replacements  = (int*) calloc (cell_count_z, sizeof (int));

        if (target_implanted_ions == NULL) return -3002;  /* cannot allocate memory */
        if (target_replacing_ions == NULL) return -3003;  /* cannot allocate memory */
        if (target_total_displacements == NULL) return -3004;  /* cannot allocate memory */
        if (target_total_interstitials == NULL) return -3005;  /* cannot allocate memory */
        if (target_total_vacancies     == NULL) return -3006;  /* cannot allocate memory */
        if (target_total_replacements  == NULL) return -3007;  /* cannot allocate memory */

        if (target_depth_implanted_ions == NULL) return -3008;  /* cannot allocate memory */
        if (target_depth_replacing_ions == NULL) return -3009;  /* cannot allocate memory */
        if (target_depth_total_displacements == NULL) return -3010;  /* cannot allocate memory */
        if (target_depth_total_interstitials == NULL) return -3011;  /* cannot allocate memory */
        if (target_depth_total_vacancies     == NULL) return -3012;  /* cannot allocate memory */
        if (target_depth_total_replacements  == NULL) return -3013;  /* cannot allocate memory */

        if (detailed_sputtering == 1) {  /* for detailed calucation of sputtering, we need these */
            target_total_sputtered  = (int*) calloc (cell_count, sizeof(int));
            if (target_total_sputtered  == NULL) return -3014;   /* cannot allocate */
        }
        if (store_energy_deposit == 1) {  /* for detailed storing of deposited energy, we need these */
            target_energy_phonons   = (double*) calloc (cell_count, sizeof (double));
            target_energy_electrons = (double*) calloc (cell_count, sizeof (double));

            target_depth_energy_phonons   = (double*) calloc (cell_count_z, sizeof (double));
            target_depth_energy_electrons = (double*) calloc (cell_count_z, sizeof (double));

            if (target_energy_phonons   == NULL) return -3015;  /* cannot allocate */
            if (target_energy_electrons == NULL) return -3016;  /* cannot allocate */

            if (target_depth_energy_phonons   == NULL) return -3017;  /* cannot allocate */
            if (target_depth_energy_electrons == NULL) return -3018;  /* cannot allocate */
        }

        /* go through materials, init arrays for interstitials and vacancies */
        for (i=0; i<number_of_materials; i++) {
            /* create arrays of pointers to arrays */
            list_of_materials[i].target_implanted_recoils_int  =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
            list_of_materials[i].target_implanted_recoils_repl =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
            list_of_materials[i].target_elemental_vacancies   =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
            list_of_materials[i].target_elemental_disp        =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);

            list_of_materials[i].target_depth_implanted_recoils_int  =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
            list_of_materials[i].target_depth_implanted_recoils_repl =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
            list_of_materials[i].target_depth_elemental_vacancies   =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
            list_of_materials[i].target_depth_elemental_disp        =
                (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);

            if (list_of_materials[i].target_implanted_recoils_int  == NULL) return -3019;
            if (list_of_materials[i].target_implanted_recoils_repl == NULL) return -3020;
            if (list_of_materials[i].target_elemental_vacancies    == NULL) return -3021;
            if (list_of_materials[i].target_elemental_disp         == NULL) return -3022;

            if (list_of_materials[i].target_depth_implanted_recoils_int  == NULL) return -3023;
            if (list_of_materials[i].target_depth_implanted_recoils_repl == NULL) return -3024;
            if (list_of_materials[i].target_depth_elemental_vacancies    == NULL) return -3025;
            if (list_of_materials[i].target_depth_elemental_disp         == NULL) return -3026;

            /* for detailed calucation of sputtering, we need these */
            if (detailed_sputtering == 1) {
                list_of_materials[i].target_sputtered_atoms =
                    (int**) malloc (sizeof (int*) * list_of_materials[i].element_count);
                if (list_of_materials[i].target_sputtered_atoms   == NULL) return -3027;
            }

            /* go through elements, make arrays for each one */
            for (j=0; j<list_of_materials[i].element_count; j++) {
                list_of_materials[i].target_implanted_recoils_int[j]  =
                    (int*) calloc (cell_count, sizeof (int));
                list_of_materials[i].target_implanted_recoils_repl[j] =
                    (int*) calloc (cell_count, sizeof (int));
                list_of_materials[i].target_elemental_vacancies[j]    =
                    (int*) calloc (cell_count, sizeof (int));
                list_of_materials[i].target_elemental_disp[j]         =
                    (int*) calloc (cell_count, sizeof (int));

                list_of_materials[i].target_depth_implanted_recoils_int[j]  =
                    (int*) calloc (cell_count_z, sizeof (int));
                list_of_materials[i].target_depth_implanted_recoils_repl[j] =
                    (int*) calloc (cell_count_z, sizeof (int));
                list_of_materials[i].target_depth_elemental_vacancies[j]    =
                    (int*) calloc (cell_count_z, sizeof (int));
                list_of_materials[i].target_depth_elemental_disp[j]         =
                    (int*) calloc (cell_count_z, sizeof (int));

                if (list_of_materials[i].target_implanted_recoils_int[j]  == NULL) return -3028;
                if (list_of_materials[i].target_implanted_recoils_repl[j] == NULL) return -3029;
                if (list_of_materials[i].target_elemental_vacancies[j]    == NULL) return -3030;
                if (list_of_materials[i].target_elemental_disp[j]         == NULL) return -3031;

                if (list_of_materials[i].target_depth_implanted_recoils_int[j]  == NULL) return -3032;
                if (list_of_materials[i].target_depth_implanted_recoils_repl[j] == NULL) return -3033;
                if (list_of_materials[i].target_depth_elemental_vacancies[j]    == NULL) return -3034;
                if (list_of_materials[i].target_depth_elemental_disp[j]         == NULL) return -3035;

                /* For detailed calucation of sputtering, we need these */
                if (detailed_sputtering == 1) {
                    list_of_materials[i].target_sputtered_atoms[j] =
                        (int*) calloc (cell_count, sizeof (int));
                    if (list_of_materials[i].target_sputtered_atoms[j]   == NULL) return -3036;
	              }
            }
        }

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=1) printf ("Memory for target arrays has been allocated.\n\n");
        /* MPI=============================================== */
#else
        if (print_level >= 1) printf ("Memory for target arrays has been allocated.\n\n");
#endif

        /* switch : 0 - CSG, 1 - FETM. */
        /* now read in the target composition file */
        /* target composition file for non-dyanmic version based on materials */
        switch (geometry_type) {
        case 0 :  /* multi_layer bulk geometry */
            result = read_bulk_shape (TargetCompositionFileName);  /* (1) - (8) - (2) */
            if (result != 0) {
                printf ("Error: cannot read multi-layer bulk geometry file.\n");
                return result;
            }
            else {
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node == ROOT) printf ("This is the multi-layer bulk geometry version of iran3d.\n\n");
                /* MPI=============================================== */
#else
                printf ("This is the multi-layer bulk geometry version of iran3d.\n\n");  // %s, shape
#endif
            }
            break;
        case 1 :  /* csg geometry */
            result = read_csg_shape (TargetCompositionFileName);  /* (1) - (8) - (2) */
            if (result != 0) {
                printf ("Error: cannot read csg geometry file.\n");
                return result;
            }
            else {
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node == ROOT) printf ("This is the CSG geometry version of iran3d.\n\n");
                /* MPI=============================================== */
#else
                printf ("This is the CSG geometry version of iran3d.\n\n");  // %s, shape
#endif
            }
            break;
        case 2 :  /* fetm geometry */
            result = read_fetm_shape (TargetCompositionFileName);  /* (1) - (8) - (3) */
            if (result != 0) {
                printf ("Error: cannot read fetm geometry file.\n");
                return result;
            }
            else {
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node == ROOT) printf ("This is the FETM geometry version of iran3d.\n\n");
                /* MPI=============================================== */
#else
                printf ("This is the FETM geometry version of iran3d.\n\n");
#endif
            }
            break;
        default :
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node == ROOT) printf ("Unknown geometry 'geometry_type':, %i", geometry_type);
            /* MPI=============================================== */
#else
            printf ("Unknown geometry 'geometry_type':, %i", geometry_type);
#endif
            break;
        }
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0)
            printf ("Target composition read from %s.\n", TargetCompositionFileName);
        /* MPI=============================================== */
#else
        if (print_level >= 0)
            printf ("Target composition read from %s.\n", TargetCompositionFileName);
#endif
    }
    else {  /* only calc memory usage */
#ifndef MPI_PRALLEL
        /* MPI=============================================== */
        lui_temp += 7 * sizeof (int) * cell_count;
        lui_temp += sizeof (float) * cell_count;
        if (detailed_sputtering == 1) {  /* for detailed calucation of sputtering, we need these */
            lui_temp += sizeof (int) * cell_count;
            lui_temp += sizeof (char) * cell_count;
        }
        /* for detailed storing of deposited energy */
        if (store_energy_deposit == 1) lui_temp += 2 * sizeof (double) * cell_count;
        for (i=0; i<number_of_materials; i++) {
            for (j=0; j<list_of_materials[i].element_count; j++) {  /* go through elements,
                                                             make arrays for each one */
                lui_temp += 4 * cell_count * sizeof (int);
                /* for detailed calucation of sputtering, we need these */
                if (detailed_sputtering == 1) lui_temp += cell_count * sizeof (int);
            }
        }
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Target historgrams:          %li bytes\n", lui_temp);
        /* MPI=============================================== */
#endif
    }

    return 0;
}

/*==== now come some helping functions concerning file reading ====*/

/*=============================================================================
  Function Name : read_target_structure_data_block
  Description   : Needs to be called from the ini file reader while the traget
                  structure input file is read.

  Inputs  :
            char* par_name
            char* par_value
  Outputs : no.

  Notes :
      We won't need this at the moment...
=============================================================================*/
int read_target_structure_data_block (char *block_name) {
    return 0;
}

/*=============================================================================
  Function Name : read_target_structure_data
  Description   : Needs to be called from the ini file reader while the target
                  structure input file is read.

  Inputs  :
            char* par_name
            char* par_value
  Outputs : no.

  Notes : no
=============================================================================*/
int read_target_structure_data (char *par_name, char *par_value) {

    if (strcmp (par_name, "cell_count_x") == 0)  /* ... */
        sscanf (par_value, "%i", &cell_count_x);
    if (strcmp (par_name, "cell_count_y") == 0)  /* ... */
        sscanf (par_value, "%i", &cell_count_y);
    if (strcmp (par_name, "cell_count_z") == 0)  /* ... */
        sscanf (par_value, "%i", &cell_count_z);

    if (strcmp (par_name, "cell_size_x")  == 0)  /* ... */
        sscanf (par_value, "%f", &cell_size_x);
    if (strcmp (par_name, "cell_size_y")  == 0)  /* ... */
        sscanf (par_value, "%f", &cell_size_y);
    if (strcmp (par_name, "cell_size_z")  == 0)  /* ... */
        sscanf (par_value, "%f", &cell_size_z);

    if (strcmp (par_name, "sub_surf_z")   == 0)  /* ..., ygli */
        sscanf (par_value, "%f", &sub_surf_z);

    if (strcmp (par_name, "CompositionFileType") == 0)  /* ... */
        sscanf (par_value, "%i", &target_composition_file_type);

    if (strcmp (par_name, "CompositionFileName") == 0)  /* ... */
        if (single_input_file != 1) strcpy (TargetCompositionFileName, par_value);

    return 0;
}

/*=============================================================================
  Function Name : target_index
  Description   : Calcualtes index for accessing one-dimensional target arrays
                  knowing the integer cell coordinates.

  Inputs  : int x, int y, int z
  Outputs : no.

  Notes :
         layer_count_xy = cell_count_x * cell_count_y,
         still use layer_count_yz is in order to make output in style of
         (x, y, z, value).
=============================================================================*/
int get_target_index (int x, int y, int z) {
    return x * layer_count_yz + y * cell_count_z + z;
}

/*=============================================================================
  TO DO: Idea: use reciprocal values for 1/cell_count and so on,
  to speed up calculations.
=============================================================================*/

/*=============================================================================
  Function Name : get_target_XYZ
  Description   : Calculates the x,y,z coordinates (cell numbers) for a given
                  index ...; so the inverse function of TargetIndex().

  Inputs  :
            int index
  Outputs : int* x, int* y, int* z

  Notes :
         layer_count_xy = cell_count_x * cell_count_y,
         still use layer_count_yz is in order to make output in style of
         (x, y, z, value).
=============================================================================*/
void get_target_XYZ (int index, int *x, int *y, int *z) {
    *z = ( index % layer_count_yz) % cell_count_z;
    *y = ((index % layer_count_yz) - (*z)) / cell_count_z;
    *x = ((index - (*z) - cell_count_z * (*y)) / layer_count_yz);

    return;
}

/*=============================================================================
  Function Name : get_relative_XYZ
  Description   : Calculates the relative x,y,z coordinates (unit: traget_size)
                  for a given position.

  Inputs  : int x, int y, int z

  Outputs : float *rx, float *ry, float *rz

  Notes : no.
=============================================================================*/
void get_relative_XYZ (int x, int y, int z, float *rx, float *ry, float *rz) {

    *rx = (float) x / (float) cell_count_x;
    *ry = (float) y / (float) cell_count_y;
    *rz = (float) z / (float) cell_count_z;

    return;
}

/*=============================================================================
  Function Name : get_cell_index
  Description   : Calculates index for accessing one-dimensional target arrays
                  from a float-value point in space, if in the target region.

  Inputs  : float x, float y, float z
  Outputs : int cell_i

  Notes : no.
=============================================================================*/
int get_cell_index (double x, double y, double z, int *cell_i) {
    int ix, iy, iz;

    if (x >= 0.0 && x < target_size_x && y >= 0.0 && y < target_size_y && z >= 0.0 &&
        z < target_size_z) {
        ix = x / cell_size_x;
        iy = y / cell_size_y;
        iz = z / cell_size_z;

        *cell_i = get_target_index (ix, iy, iz);  /* Determine initial cell */
        return 1;
    }
    else {
        return 0;
    }
}

/*=============================================================================
  Function Name : cal_relative_target_atom_position
  Description   : This calculates the direction in which the target nucleus is found.
                  v is the projectile velocity vector, the IP vector is returned
                  in p components.

  Inputs  :
            float vx,float vy, float vz
            unsigned int azim_angle
  Outputs : float *px, float *py, float *pz

  Notes :
      This routine works similar to the rotation (DIRCOS) routine. It simply
      assumes a fixed scattering angle of 90 degrees and reverses the resulting
      vector. Adding the result multiplied with impact_par to the current
      projectile position leads to the target nucleus position.

      Function call:
        fromcorteo - float my_inv_sqrt (float val);
                     float sqrtdf (double val).
=============================================================================*/
void cal_relative_target_atom_position (double vx, double vy, double vz, double *px,
                             double *py, double *pz, unsigned int iazim_angle) {
    float k, kinv;
    float k2 = 1.0f - vz * vz;

    /* random azimutal rotation angle components */
    float cosomega = cos_azim_angle[iazim_angle];
    float sinomega = sin_azim_angle[iazim_angle];
    /* In the rotation routine, we would have to increase the azim-counter here, but
       since we want to have the same angle for target position and deflection, we must
       not change the index in this function! */
    if (k2 < 0.000001) {  /* extremely rare case */
        *pz = 0;
        *py = cosomega;
        *px = sinomega;
        return;
    }

    /* usual case */
#ifndef SAFE_SQR_ROTATION
    kinv = inv_sqrt (k2);  /* 1/sqrt() (approximate) */
    k = k2 * kinv;         /* so we get k and 1/k by a table lookup and a multiplication */
#else
    /* these two lines can be replaced by the following two lines but... */
    k  = sqrtdf (k2);   /* ... using a sqrt() here makes the program 25% slower! */
    kinv = 1.0f / k;
#endif

    *px = -kinv * (vx * vz * cosomega + vy * sinomega);
    *py = -kinv * (vy * vz * cosomega - vx * sinomega);
    *pz =  k * cosomega;

#ifdef SAFE_ROTATION  /* makes slower, but safer */
    if (*px > 1)  *px =  1;
    if (*px < -1) *px = -1;
    if (*py > 1)  *py =  1;
    if (*py < -1) *py = -1;
    if (*pz > 1)  *pz =  1;
    if (*pz < -1) *pz = -1;
#endif

    return;
}
