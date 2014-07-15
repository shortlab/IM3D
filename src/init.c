/******************************************************************************
  Module Name : init.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Contains some initialization functions.

  Others :
      Error numbers in this file: 1000 - 1999.
      Refer to iradina and Corteo.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "init.h"

/* define the type screening function, can be moved into Config_in */
#define screening_type 1  /* Screening function type used to generate corteo.*.mat:
                             0 - unscreened potential (screening = 1);
                             1 - universal screening function;
                             2 - Lenz-Jensen screening function;
                             others - unknown. */

/*=============================================================================
  Function Name : init_configuration
  Description   : Read configuration from file, initialize variables etc.

  Inputs  : char* ConfigFileName
  Outputs : no.

  Notes :
      Function call:
        fileio - int check_split_input_file (char* filename);
                 int read_init_file (int(*DataBlockReader) (char* BlockName),
                                     int (*DataReader) (char* ParName,
                                     char* ParValue), char* filename);
                 int load_inverse_Erf ();
        fromcorteo - int load_matrix (void);
                     void compute_lists ();
                     void fill_my_sqrt_table ();
                     void randomize_list (float *list, unsigned int maxlist);
                     float d2f (double val);
                     double randomx ();
        utils - int load_Chu_straggling_values ();
                int calculate_normalization_factor (int num_of_ions);
                int calculate_normalization_factor (int num_of_ions);
        target - int initialize_materials (char* Filename);
                 int init_target_structure (char* Filename).
=============================================================================*/
int init_configuration (char *ConfigFileName) {
    int   i, result, exist;
    float length, random;
#ifndef MPI_PRALLEL
    /* MPI=============================================== */
    long unsigned int lui_temp;
    /* MPI=============================================== */
#endif
    struct transmitted_ion temp;
    char file_name[100];

    /* some default values */
    simulation_type  = 0;
    straggling_model = 0;

    /* if the config file is a combined input file (including target, materials,
       and composition), then it must be split up into the separate files first: */
    check_split_input_file (ConfigFileName);

    /* read general config: */  /* (1) - (1) */
    result = read_init_file (read_config_file_data_block, read_config_file_data, ConfigFileName);
    if (result != 0) {
        printf("Error reading config file %s.\n", ConfigFileName);
        return result;
    }
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    /* MPI random */
    seed1 = mpi_seed (seed1);
    seed2 = mpi_seed (seed2);
    printf ("Seed in node %i is: %i\t%i\n", my_node, seed1, seed2);

    if (my_node==ROOT && print_level>=0)
        printf ("Configuration read from %s.\n", ConfigFileName);
    /* MPI=============================================== */
#else
    if (print_level >= 0) printf ("Configuration read from %s.\n", ConfigFileName);
#endif

    if (mem_usage_only == 0) {  /* do real stuff, not just estimating memory usage */
        /* calculate or load the corteo scattering matrix */
        switch (screening_type) {
        case 0 :  /* unscreened potential (screening = 1) */
            sprintf(file_name, "data/corteo.none.mat");
            break;
        case 1 :  /* universal screening function */
            sprintf(file_name, "data/corteo.univ.mat");
            break;
        case 2 :  /* Lenz-Jensen screening function */
            sprintf(file_name, "data/corteo.none.mat");
            break;
        default :  /* unknown */
            fprintf(stderr, "Error: screening function of type %u unknown. Stopping here.",
                    screening_type);
            exit(1);
        }

        exist = access (file_name, 0);
        /* if file is not exist, calculate it */
        if (exist == -1) calc_matrix (screening_type, 1, file_name);
        result = load_matrix (file_name);  /* (1) - (2) */
        if (result != 1) {printf ("Error load corteo scattering matrix.\n"); return -1001;}
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf ("Corteo scattering matrix loaded.\n");
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf ("Corteo scattering matrix loaded.\n");
#endif

        /* call corteo's list generator */
        compute_lists ();  /* (1) - (3) */
        fill_fast_sqrt_table ();  /* (1) - (4) */
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf ("Lists of random numbers generated.\n");
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf ("Lists of random numbers generated.\n");
#endif

        /* normalize ion velocity unit vector */
        length = 1.0f / sqrt (ion_vx * ion_vx + ion_vy * ion_vy + ion_vz * ion_vz);
        ion_vx *= length;
        ion_vy *= length;
        ion_vz *= length;

        /* make sure that everything is correctly initialzed: */
        if (OutputFileBaseName == NULL) {
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node == ROOT)
                printf ("No output file basename specified. Using: default_out/out.\n");
            /* MPI=============================================== */
#else
            printf ("No output file basename specified. Using: default_out/out.\n");
#endif
            OutputFileBaseName = (char*) malloc (sizeof (char) * 18);
            strcpy (OutputFileBaseName, "default_out/out");
        }

        /* init array for storing transmitted ions */
        if (store_transmitted_ions == 1) {
            transmission_pointer = 0;
            transmit_list = malloc (sizeof (temp) * (max_no_ions + 1));
            if (transmit_list == NULL) {
                printf ("Not enough memory to store transmitted ions!\n");
                return -1002;
            }
        }

        /* read Chu's straggling data */
        result = load_Chu_straggling_values ();  /* (1) - (5) */
        if (result != 0) {printf ("Error reading Chu's straggling values!\n"); return result;}
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf ("Chu's straggling data read.\n");
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf ("Chu's straggling data read.\n");
#endif

        /* read invserse error function list and randomize it*/
        result = load_inverse_Erf ();  /* (1) - (6) */
        if (result != 0) { printf ("Error reading invser Erf() list!\n"); return result;}
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf ("Invsere Erf list read.\n");
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf ("Invsere Erf list read.\n");
#endif
        randomize_list (inverse_erf_list, MAXERFLIST);
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf ("Invsere Erf list randomized.\n");
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf ("Invsere Erf list randomized.\n");
#endif
    }
    else {  /* Esimate memory usage: */
#ifndef MPI_PRALLEL
        /* MPI=============================================== */
        lui_temp = (DIME * DIMS * sizeof (float));
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Scattering matrix (reduced): %li bytes\n", lui_temp);
        lui_temp = MAXRANLIST * sizeof (float) * 2 + MAXLOGLIST * sizeof (float) * 2
                 + MAXAZILIST * sizeof (float) * 2;
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Random lists:                %li bytes\n", lui_temp);
        lui_temp = 65536 * 2;  /* sqrt lists */
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Sqrt lists:                  %li bytes\n", lui_temp);
        lui_temp = sizeof (temp) * (max_no_ions + 1);
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Transmitted ions:            %li bytes\n", lui_temp);
        lui_temp = 400 * sizeof (float);
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Chu straggling:              %li bytes\n", lui_temp);
        lui_temp = MAXERFLIST * sizeof (float);
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY ERF list:                    %li bytes\n", lui_temp);
#endif
    }

    /* read and init materials or element data: */
    result = init_materials (MaterialsFileName);  /* (1) - (7) */
    if (result != 0) {
        printf ("Error reading materials file %s.\n", MaterialsFileName);
        return result;
    }
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=0)
        printf ("Materials read from %s.\n", MaterialsFileName);
    /* MPI=============================================== */
#else
    if (print_level >= 0) printf ("Materials read from %s.\n", MaterialsFileName);
#endif
    for (i=0; i<8; i++) {  /* empty ion leaving counter */
        leaving_ions[i] = 0;
    }

    /* read and init target structure */
    result = init_target_structure (TargetStructureFileName);  /* (1) - (8) */
    if (result != 0) {
        printf ("Error reading target structure file %s.\n", TargetStructureFileName);
        return result;
    }
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=0)
        printf ("Target structure read from %s.\n", TargetStructureFileName);
    /* MPI=============================================== */
#else
    if (print_level >= 0)
        printf ("Target structure read from %s.\n", TargetStructureFileName);
#endif

    /* init some random values (random starting point for lists of random numbers */
    random = d2f (randomx ());
    azim_angle       = (unsigned int) (random * MAXAZILIST);
    ran_list         = (unsigned int) (random * MAXRANLIST);
    ran_log_list     = (unsigned int) (random * MAXLOGLIST);
    erf_list_pointer = (unsigned int) (random * MAXLOGLIST);

    /* init array for ion single ion sputter yield histogram */
    if (single_ion_sputter_yields == 1) {
        sputter_yield_histogram = malloc ((MAX_SPUTTERED + 1) * sizeof (int));
        if (sputter_yield_histogram == NULL) {
            printf ("Error: Not enough memory store sputter histogram!\n");
            return -1003;
        }
        for (i=0; i<MAX_SPUTTERED; i++) sputter_yield_histogram[i] = 0;
    }

    calculate_normalization_factor (max_no_ions);
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=0)
        printf ("Normalization factor:\t%lg\n", unit_conversion_factor);
    /* MPI=============================================== */
#else
    if (print_level >= 0) printf ("Normalization factor:\t%lg\n", unit_conversion_factor);
#endif

    return 0;
}

/******************************************************************************
  Next 3 functions: table-based sqrt and 1/sqrt evaluation
  tables are computed for the sqrt and 1/sqrt of mantissa and exponent
  separately (one table, the shortest one related to the exponent, could
  be eliminated by simply bit-shifting the exponent by one bit... to be explored).
******************************************************************************/
/*=============================================================================
  Function Name : fill_fast_sqrt_table
  Description   : Fill tables.

  Inputs  : no.
  Outputs : no.

  Notes :
      Function call:
        fromcorteo - float sqrtdf (double val).
=============================================================================*/
void fill_fast_sqrt_table (void) {
    unsigned long i, j, n = 1<<16;  /* 1<<16: mantissa tables contain 65536 values */
    float val;

    /* mantissa from 0 to 2 in 65536 steps, giving off precision! */
    for (i=0; i<n; i++) {
        /* generate a float values */
        j = (i<<7) + ((1<<23) * 127);
        val = *(float *)&j;

        /* store its sqrt and 1/sqrt in tables */
        sqrt_table[i] = sqrtdf (val);
        inv_sqrt_table[i] = 1.0f / sqrt_table[i];
    }

    /* exponent from 2^-128 to 2^128 */
    for (i=0; i<(1<<8); i++) {
        /* generate a float values */
        j = i<<23;
        val = *(float *)&j;

        /* store its sqrt and 1/sqrt in tables */
        sqrt_table_exp[i]    = sqrtdf (val);
        inv_sqrt_table_exp[i] = 1.0f / sqrt_table_exp[i];
    }

    return;
}

/*=============================================================================
  Function Name : load_Chu_straggling_values
  Description   : Loads the Chu straggling data; returns 0 on success.
                  (from Q.Yang et al., NIMB vol. 61, page 149, (1991)).

  Inputs  : no.
  Outputs : no.

  Notes :
      This function is adapted from corteo.

      Function call:
        utils - void ignore_line (FILE *ifp).
=============================================================================*/
int load_Chu_straggling_values () {
    FILE *fp;
    unsigned int k, l, Z;
    char temp[1000];

    fp = fopen ("data/chu.dat", "r");
    if (fp == NULL) return -1004;

    ignore_line (fp);  /* skip first line */
    for (k=0; k<92; k++) {
        fscanf (fp, "%u", &Z);
        for (l=0; l<4; l++) {
            fscanf (fp, "%s", temp);
            chu_values[Z][l] = a2f (temp);
        }
        if (Z != k+2) return -1005;  /* check consistency */
    }
    fclose(fp);

    return 0;
}

/*=============================================================================
  Function Name : load_inverse_Erf
  Description   : Load a list of inverse error function erfinv(x) values that
                  contains MAXERFLIST elements the list is such that x is
                  uniformly distributed between -1 and 1 excluding these
                  boundaries i.e. x = -1+dx, -1+2dx, ...1-2dx, 1-dx, with
                  1/2dx = MAXERFLIST.
                  Return 0 on success, some other value else.

  Inputs  : no.
  Outputs : no.

  Notes :
      Function adapted from corteo.c.

      function call:
        fromcorteo - float a2f (char * s).
=============================================================================*/
int load_inverse_Erf () {
    unsigned int k;
    char erf_val[1000];
    FILE *ifp;

    ifp = fopen ("data/erfinv.dat", "r");
    if (ifp == NULL) return -1006;

    for (k=0; k<MAXERFLIST; k++) {
        fscanf (ifp, "%s", erf_val);
        inverse_erf_list[k] = a2f (erf_val);
    }
    fscanf (ifp, "%s", erf_val);
    fclose (ifp);

    /* control value at the end should be the number of elements in the list */
    if (atoi (erf_val) != MAXERFLIST) {
        return -1007;
    }

    return 0;
}
