/******************************************************************************
  Module Name : im3d.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Main program.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
      12/06/13            started designing;
      02/26/14            started coding;
      04/29/14            finished version_1.0.0;
      11/10/14            finished version_1.0.4;
      11/26/14            renamed as im3d_static_v1.0;
      01/08/15            Revise some subroutines;
      02/27/15            Add the radial distribution functions statistics;
      05/12/15            version im3d_static_v1.1.0;
      05/14/15            Add the modulus aivxyz;
******************************************************************************/
#include "im3d.h"

/*=============================================================================
  Function Name : main
  Description   : Main program.

  Inputs  :
            int argc - number of commonds;
            char* argv[] - commonds.
  Outputs : no.

  Notes :
      Function call:
        im3d - int display_startup_message ();
        utils - int handle_cmd_line_options (int argc, char* argv[]);
                int write_status_file(char* status_text, int ion_number);
                int init_configuration (char* ConfigFileName);
                int calculate_normalization_factor (int num_of_ions);
                int store_results (char* BaseName);
                int convert_material_to_element (char* OutputFile);
        transport - int irradiate_target ().
=============================================================================*/
int main (int argc, char *argv[]) {
    int result;
    int cmd_result;
    clock_t start_time = 0, finish_time;
    float   duration_time;

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_node  );
    MPI_Comm_size (MPI_COMM_WORLD, &num_nodes);
    if (my_node == ROOT) start_time = MPI_Wtime ();
    /* MPI=============================================== */
#else
    /* run time start */
    start_time = clock ();
#endif

    /* init some default values for some global variables */
    print_level              = 0;
    wait_before_end          = 0;
    mem_usage_only           = 0;
    mem_usage_details        = 0;
    mem_usage                = 0;
    override_max_ions        = 0;
    override_energy          = 0;
    display_interval         = 100;
    store_path_limit         = -1;
    store_exiting_recoils    = 0;
    store_exiting_limit      = 10;
    tracing_recoil_or_not    = 1;  /* KP model, Aug. 05, 2014 */
    scattering_calculation   = 0;
    max_annular_coll_volumes = 0;
    transport_type           = 0;
    normalize_output         = 0;  /* default: no normalization */
    do_not_store_damage      = 0;
    unit_conversion_factor   = 1.0;
    create_status_file       = 0;
    single_input_file        = 0;
    hydrogen_in_target       = 0;
    status_update_interval   = 100000;
    conv_create_separate_elements = 0;

    /* obtain memory for file_names */
    ConfigFileName = (char*) malloc (sizeof (char) * 1024);
    if (ConfigFileName == NULL) {printf ("Not enough memory\n"); return -1;}
    strcpy (ConfigFileName, "Config.in");  /* Default */
    MaterialsFileName = (char*) malloc (sizeof (char) * 1024);
    if (MaterialsFileName == NULL) {printf ("Not enough memory\n"); return -1;}
    strcpy (MaterialsFileName, "Materials.in");  /* Default */
    ConversionFileName = (char*) malloc (sizeof (char) * 1024);
    if (ConversionFileName == NULL) {printf ("Not enough memory\n"); return -1;}
    strcpy (ConversionFileName, "Converted.input");  /* Default */
    ElementsFileName = (char*) malloc (sizeof (char) * 1024);
    if (ElementsFileName == NULL) {printf ("Not enough memory\n"); return -1;}
    strcpy (ElementsFileName, "Elements.in");  /* Default */

    TargetStructureFileName = (char*) malloc (sizeof (char) * 1024);
    if (TargetStructureFileName == NULL) {printf ("Not enough memory\n"); return -1;}
    TargetCompositionFileName = (char*) malloc (sizeof (char) * 1024);
    if (TargetCompositionFileName == NULL) {printf ("Not enough memory\n"); return -1;}

    /* handle command line arguments: */
    cmd_result = handle_cmd_line_options (argc, argv);
    if (cmd_result == 0) {  /* normal simulation */
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node == ROOT) {
            display_startup_message ();
            if (create_status_file == 1) write_status_file ("init0", 0);
        }
        /* MPI=============================================== */
#else
        /* feedback and display startup message */
        display_startup_message ();
        if (create_status_file == 1) write_status_file ("init0", 0);
#endif

        /* init all things */
        result = init_configuration (ConfigFileName);  /* (1) */
        if (result != 0) {
            printf ("Initialization error: %i. Aborting.\n", result);
            return result;
        }
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && create_status_file==1) write_status_file ("init1", 0);
        /* MPI=============================================== */
#else
        if (create_status_file == 1) write_status_file ("init1", 0);

        /* aiv_xyz */
        init_aivxyz ();
#endif

        /* do not simulate, only show memory usage. estimate about 1 MB overhead
           (for arrays of pointers, recursive function calls and so on) */
#ifndef MPI_PRALLEL
            /* MPI=============================================== */
        if (mem_usage_only == 1) {
            printf ("Estimated total memory requirement: %li MByte\n", (mem_usage / 0x100000) + 1);
            return 99;
            /* MPI=============================================== */
        }
#endif

        /* the command line option has been used to override ion number */
        if (override_max_ions > 0) {
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node==ROOT && print_level>=1)
                printf ("Maximum number of ions limited to %i.\n", override_max_ions);
            /* MPI=============================================== */
#else
            if (print_level >= 1)
                printf ("Maximum number of ions limited to %i.\n", override_max_ions);
#endif
            max_no_ions = override_max_ions;
        }

        /* the command line option has been used to override ion energy */
        if (override_energy > 0.0000001) {
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node==ROOT && print_level>=0)
                printf ("New ion energy: %g.\n", override_energy);
            /* MPI=============================================== */
#else
            if (print_level >= 0) printf ("New ion energy: %g.\n", override_energy);
#endif
            ion_initial_energy = override_energy;
        }

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node == ROOT) {
            if (create_status_file == 1) write_status_file ("init2", 0);
            if (print_level>=1){
                printf ("\nEverything is initialized and ready.\n");
                printf ("-----------------------------------------------\n\n");
            }
            if (print_level >= -1) printf ("\nStarting simulation of irradiation...\n");
        }
#else
        if (create_status_file == 1) write_status_file ("init2", 0);
        if (print_level>=1){
            printf ("\nEverything is initialized and ready.\n");
            printf ("-----------------------------------------------\n\n");
        }
        if (print_level >= -1) printf ("\nStarting simulation of irradiation...\n");
#endif

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node == ROOT && print_level>=1) {
            calculate_normalization_factor (max_no_ions);
            printf ("Normalization factor is: %g\n", unit_conversion_factor);
            fflush (stdout);
        }
        /* MPI=============================================== */
#else
        if (print_level >= 1) {
            calculate_normalization_factor (max_no_ions);
            printf ("Normalization factor is: %g\n", unit_conversion_factor);
            fflush (stdout);
        }
#endif

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        mpi_init_array ();

        mpi_distribute ();

        MPI_Barrier (MPI_COMM_WORLD);
        /* MPI=============================================== */
#endif

        /* irradiate the target with ions */
        irradiate_target ();  /* (2) */

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        MPI_Barrier (MPI_COMM_WORLD);

        mpi_reduce_data ();
        /* MPI=============================================== */
#endif

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node == ROOT) {
            if (print_level >= 1) {
                printf ("\nSimulation is finished.\n");
                printf ("-----------------------------------------------\n");
            }
            if (print_level >= -1) {printf ("Storing final results: ...\n "); fflush (stdout);}
            if (create_status_file == 1) write_status_file ("simend", 0);
            calculate_normalization_factor (max_no_ions);

            switch (output_format) {
            case 0 :  /* iradina */
                store_results_iradina (OutputFileBaseName);  /* (3) */
                break;
            case 1 :  /* cfg */
                store_results_cfg (OutputFileBaseName);
                break;
            case 2 :  /* msh */
                store_results_msh (OutputFileBaseName);
                break;
            case 3 :  /* vtk */
                store_results_vtk (OutputFileBaseName);
                break;
            default :
                printf ("ERROR: output format is not existent.");
                break;
            }

            if (print_level >= -1) printf ("done.\n\n");
            if (wait_before_end == 1) {printf ("Press enter to exit.\n"); getc (stdin);}
            if (create_status_file == 1) write_status_file ("end", 0);
        }
        /* MPI=============================================== */
#else
	    if (print_level >= 1) {
            printf ("\nSimulation is finished.\n");
            printf ("-----------------------------------------------\n");
        }
        if (print_level >= -1) {printf ("Storing final results: ...\n "); fflush (stdout);}
        if (create_status_file == 1) write_status_file ("simend", 0);
        calculate_normalization_factor (max_no_ions);

        switch (output_format) {
        case 0 :  /* iradina */
            store_results_iradina (OutputFileBaseName);  /* (3) */
            break;
        case 1 :  /* cfg */
            store_results_cfg (OutputFileBaseName);
            break;
        case 2 :  /* msh */
            store_results_msh (OutputFileBaseName);
            break;
        case 3 :  /* vtk */
            store_results_vtk (OutputFileBaseName);
            break;
        default :
            printf ("ERROR: output format is not existent.");
            break;
        }

        /* aiv_xyz */
        store_aivxyz_cfg (OutputFileBaseName);

        if (print_level >= -1) printf ("done.\n\n");
        if (wait_before_end == 1) {printf ("Press enter to exit.\n"); getc (stdin);}
        if (create_status_file == 1) write_status_file ("end", 0);
#endif
    }
    else {  /* do something else than simulate */
#ifndef MPI_PRALLEL
        /* MPI=============================================== */
        if (cmd_result < 0) {  /* an error occurred */
            return cmd_result;
        }
        else {  /* perform some operation that is not simulation */
            switch (cmd_result) {
            case 3 :  /* convert material to element file */
                printf ("\nYou are using im3d to convert a material based target definition to an element based one.\n\n");

                /* init all things, needs to be done to read material file */
                result = init_configuration (ConfigFileName);

                if (result != 0) {
                    printf ("Initialization error: %i. Aborting.\n", result);
                    return result;
                }

                /* call converter: */
                result = convert_material_to_element (ConversionFileName);  /* (4) */
                if (result != 0) {
                    printf ("Conversion failed. Error: %i\n", result);
                    return result;
                }

                break;
            case 1 :
                break;  /* help was printed. Do nothing else */
            default :
                printf ("Error: Unknown operation to perform.\n");
                break;
            }

        }
        /* MPI=============================================== */
#endif
    }

    /* run time finish */
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node == ROOT) {
        finish_time = MPI_Wtime ();
        duration_time = (float) (finish_time - start_time);
        printf ("Run time: %f seconds.\n", duration_time);
    }

    MPI_Finalize ();
    /* MPI=============================================== */
#else
    finish_time = clock ();

    /* run time */
    duration_time = (float) (finish_time - start_time) / CLOCKS_PER_SEC;
    printf ("Run time: %f seconds.\n", duration_time);
#endif

    return 0;
}

/*=============================================================================
  Function Name : display_startup_message
  Description   : Feedback and display startup message.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
int display_startup_message () {
    if (print_level >= -1) {
        printf ("************************************************************************\n");
        printf ("IM3D: Ion Irradiation of nanostructured Materials\n");
        printf ("       -- a 3D parallel Monte Carlo simulation code\n\n");
        printf ("IM3D Version %i.%i.%i%s, ", VERSION, SUBVERSION, SUBSUBVERSION, RELEASESTRING);
        printf ("%s\n%s\n", VERSIONDATE, VERSIONCOMMENT);
        printf ("by Yonggang Li (Y.G. Li), 2014, ygli@theory.issp.ac.cn & ygli@mit.edu;\n");
        printf ("Institute of Solid Status of Physics, Chinese Academy of Sciences;\n");
        printf ("& Nuclear Science and Engineering, Massachusetts Institute of Technology.\n\n");
        printf ("Refers to:\n  H.M. Li, H.Y. Wang, Y.G. Li & Z.J. Ding*'s CSG/FETM geometry methods;\n");
        printf ("& C. Borschel's OSS 'iradina' and F. Schiettekatte's OSS 'corteo'.\n\n");
        printf ("This program comes with ABSOLUTELY NO WARRANTY.\n");
        printf ("This is a free software, and you are welcome to redistribute it under\n");
        printf ("certain conditions, see LICENCE for details.\n");
        printf ("************************************************************************\n\n");
    }

    return 0;
}
