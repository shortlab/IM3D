/******************************************************************************
  Module Name : material.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Contains the material-related functions, etc.

  Others :
      Error codes used in this modul: 2000 - 2999.
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "material.h"

/*=============================================================================
  Function Name : init_materials
  Description   : Read materials from input file and do preparatory calculations, etc.
                  It must be called before calling InitializeTarget.
                  Returns 0 on success.

  Inputs  : char* file_name
  Outputs : no.

  Notes :
      TODO: free arrays
      Function call:
        fileio - int read_init_file (int (*read_data_block) (char* block_name),
                    int (*read_data) (char* par_name, char* par_value),
                    char* file_name);
        fromcorteo - float d2f (double val);
                     float sqrtdf (double val);
        utils - int prepare_stopping_tables ();
                int prepare_straggling_tables (int model);
                int count_existing_elements (int* elementarray);
        target - int prepare_scattering_matrix (struct scattering_matrix * ScatMatrix,
                                float ProjM, int ProjZ, float TarM, int TarZ).
=============================================================================*/
int init_materials (char *file_name) {
    int   i, j = 0, k = 0;
    int   result;
#ifndef MPI_PRALLEL
    /* MPI=============================================== */
    int   number_exist_element;  /* number of existing elements */
    unsigned long int lui_temp;
    /* MPI=============================================== */
#endif
    float temp;
    struct scattering_matrix *psm;  /* pointer to a scattering matrix */

    /* init some material arrays etc. */
    number_of_materials = 0;
    for (i=0; i<MAX_NO_MATERIALS; i++)
        sprintf (list_of_materials[i].name, "new material %i", i);  /* some initial name */

    /* read in file with material info */  /* (1) - (7) - (1) */
    result = read_init_file (read_materials_data_block, read_materials_data, file_name);
    if (result != 0) return result;

    /* fill existing elements with 0 */
    for (i=0; i<MAX_ELEMENT_NO; i++) existing_elements[i] = 0;

    /* loop through materials, for each:
       - normalize elemental concentrations to 1,
       - calculate mean inter-atomic distance,
       - calculate maximum impact parameter (not reduced one),
       - determine elements existing in target.
    */
    number_of_all_target_elements = 0;
    for (i=0; i<number_of_materials; i++) {
        temp = 0;
        list_of_materials[i].cell_count = 0;
        for (k=0; k<(MAX_EL_PER_MAT*8); k++)  /* reset sputter-counters: */
            list_of_materials[i].sputter_counter[k] = 0;
        for (j=0; j<list_of_materials[i].element_count; j++) {
            /* Sum up. Also indicate that this element exists in the target. */
            temp += list_of_materials[i].elements_conc[j];
            existing_elements[list_of_materials[i].elements_Z[j]] = 1;
            number_of_all_target_elements++;
        }

        /* if ion SBE is not given, calculate mean value: */
        if (list_of_materials[i].ion_surf_energy < 0) {
            list_of_materials[i].ion_surf_energy = 0;
            for (j=0; j<list_of_materials[i].element_count; j++)  /* sum up */
                list_of_materials[i].ion_surf_energy += list_of_materials[i].elements_surf_energy[j];
            list_of_materials[i].ion_surf_energy /= list_of_materials[i].element_count;
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node==ROOT && print_level>0)
                printf ("No ion SBE for material %i defined. Trying mean value: %f eV\n",
                        i, list_of_materials[i].ion_surf_energy);
            /* MPI=============================================== */
#else
            if (print_level > 0)
                printf ("No ion SBE for material %i defined. Trying mean value: %f eV\n",
                        i, list_of_materials[i].ion_surf_energy);
#endif
        }
        if (isnan (list_of_materials[i].ion_surf_energy)) {
            list_of_materials[i].ion_surf_energy = 3;
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node == ROOT) printf ("No ion SBE for material %i defined. Using %f eV\n",
                                      i, list_of_materials[i].ion_surf_energy);
            /* MPI=============================================== */
#else
            printf ("No ion SBE for material %i defined. Using %f eV\n",
                    i, list_of_materials[i].ion_surf_energy);
#endif
        }

        list_of_materials[i].mean_Z = 0;
        list_of_materials[i].mean_M = 0;
        list_of_materials[i].mean_Ed = 0;  /* Modified Kinchin-Pease Model, Aug. 05, 2014 */

        /* go through elements, normalize concentration, check for hydrogen existence*/
        for (j=0; j<list_of_materials[i].element_count; j++) {
            list_of_materials[i].elements_conc[j] = list_of_materials[i].elements_conc[j] / temp;
            /* calculate average M,Z */
            list_of_materials[i].mean_Z += list_of_materials[i].elements_conc[j]
                                         * list_of_materials[i].elements_Z[j];
            list_of_materials[i].mean_M += list_of_materials[i].elements_conc[j]
                                         * list_of_materials[i].elements_M[j];
            /* Modified Kinchin-Pease Model, Aug. 05, 2014 */
            list_of_materials[i].mean_Ed += list_of_materials[i].elements_conc[j]
                                          * list_of_materials[i].elements_disp_energy[j];

            /* Hydrogen is really within target! */
            if (list_of_materials[i].elements_Z[j] == 1) hydrogen_in_target = 1;
            /* An element like the ion is within target! */
            if (list_of_materials[i].elements_Z[j] == ion_Z) ionZ_in_target = 1;
        }

        /* For testing purposes and comparisons for SRIM, we calculate the mean screening
           length and energy reduction factor (as in ZBL85). */
        list_of_materials[i].mean_A = d2f (0.1f * SCREENCONST
            / (pow (ion_Z, 0.23) + pow (list_of_materials[i].mean_Z, 0.23)) );
        list_of_materials[i].mean_F = 10.0f * list_of_materials[i].mean_A * list_of_materials[i].mean_M
            / (ion_Z * list_of_materials[i].mean_Z * (ion_M + list_of_materials[i].mean_M) * E2);

        list_of_materials[i].mean_min_red_transfer = min_energy * list_of_materials[i].mean_F
            * (ion_M + list_of_materials[i].mean_M) * (ion_M + list_of_materials[i].mean_M)
            / (4 * ion_M * list_of_materials[i].mean_M) ;

        list_of_materials[i].density_NM          =
            list_of_materials[i].density * 1e-21;  /* calculate density also in at/nm^3 */
        list_of_materials[i].atomic_distance     =  /* for conversion from cm to nm */
            1.0 / pow (4.0 * PI * list_of_materials[i].density / 3.0, 1.0 / 3.0) * 1e7;
        list_of_materials[i].rev_atomic_distance  = 1.0 / list_of_materials[i].atomic_distance;  /* ygli */
        list_of_materials[i].layer_distance      =  /* for conversion from cm to nm */
            1.0 / pow (list_of_materials[i].density, 1.0 / 3.0) * 1e7;
        //list_of_materials[i].sqrt_atomic_distance = sqrtdf (list_of_materials[i].atomic_distance);  //?
        list_of_materials[i].sqrt_rec_fl_density   =
            1.0 / sqrt(PI * flight_length_const * list_of_materials[i].density_NM);

        list_of_materials[i].mean_impact_par = 1.0f /
            sqrtdf (PI * list_of_materials[i].density_NM * list_of_materials[i].atomic_distance);
        /* Later, we need to multiply this by:
           1. sqrt of a random number between 0 and 1   "r2" in corteo manual, select impact par randomly
           2. 1/ln(-random)                             "r1" selects flight length according to Poisson statistics
           3. 1/screening_length                         to get reduced quantity
           in order to get the actual reduced impact parameter */
    }

    existing_elements[ion_Z] = 1;  /* the ion occurs as a possible projectile of course */
    /*TO DO: Can we reduce this? For the ion we only need this for stopping... not for the ScatMatrix*/
    existing_elements[1]    = 1;  /* hydrogen stopping is needed to calculate charge state */
    /*TO DO: keep different tables of existing elements for stopping and scattering, to save memory*/

    if (mem_usage_only == 0) {
        /* Load stopping tables for electronic stopping for all possible elements in all materials */
        result = prepare_stopping_tables ();  /* (1) - (7) - (2) */
        if (result !=0) {
            switch (result) {
            case -1 :
                printf ("Error: Not enough memory for stopping tables! (%i)\n", result);
                break;
            case -2 :
                printf ("Error: Cannot read data from stopping file (%i)!\n", result);
                break;
            case -3 :
                printf ("Error: control data mismatch in stopping file %i\n", j);
                break;
            default :
                printf ("Error: Cannot load stopping tables! (%i)\n", result);
            }
            return result;
        }

        /* Load straggling tables for electronic energy loss straggling for all
           possible elements in all materials. */
        result = prepare_straggling_tables (straggling_model);  /* (1) - (7) - (3) */
        if (result != 0) {
            printf ("Error: Cannot create straggling tables!\n");
            return result;
        }

        /* Modified Kinchin-Pease model, Aug. 5, 2014 */
        if (tracing_recoil_or_not == 0) {
            //result = prepare_KP_tables1 ();
            result = prepare_KP_tables2 ();
            if (result != 0) {
                printf ("Error: Cannot create Kinchin-Pease tables!\n");
                return result;
            }
        }

        /* Go through possible target element and create scattering matrices for
	       scattering of ion from each target element. */
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) {
            printf ("Prepare scattering matrices for ion on target collisions...");
            fflush (stdout);
        }
        /* MPI=============================================== */
#else
        if (print_level >= 0) {
            printf ("Prepare scattering matrices for ion on target collisions...");
            fflush (stdout);
        }
#endif
        for (i=0; i<MAX_ELEMENT_NO; i++) {
            if (existing_elements[i] == 1) {  /* possible target element */
                psm = &(ion_scattering_matrix[i]);
                result = prepare_scattering_matrix (psm, ion_M, ion_Z, most_abundant_isotope[i], i);  /* (1) - (7) - (4) */
                /* Note! This scattering is only valid for the most abundant isotope!
                   If scattering for different isotopes was accounted this would be
                   a bit heavy on memory.
                   If more accuracy is need, we should use MAGIC instead for each collision event */
                if (result != 0) {
                    printf ("Error: insufficient memory for ion-target scattering matrices! %i \n",
                            result);
                    return result;
                }
            }
        }
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf (" finished\n"); fflush (stdout);
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf (" finished\n"); fflush (stdout);
#endif

        /* Go through all possible other projectiles (which are the target elements) and create
           scattering matrices for scattering of the projectile with any other target element. */
        /* Beware, this two dimensional matrix might get rather large... */
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) {
            printf ("Prepare scattering matrices for recoil on target collisions...");
            fflush (stdout);
        }
        /* MPI=============================================== */
#else
        if (print_level >= 0) {
            printf ("Prepare scattering matrices for recoil on target collisions...");
	        fflush (stdout);
        }
#endif
        for (j=0; j<MAX_ELEMENT_NO; j++) {
            if (existing_elements[j] == 1) {  /* possible projectile element */
                for (i=0; i<MAX_ELEMENT_NO; i++) {
                    if (existing_elements[i] == 1) {  /* possible target element */
                        result = prepare_scattering_matrix (&(scattering_matrices[j][i]),
                                 most_abundant_isotope[j], j, most_abundant_isotope[i], i);
                        /* Not: First index in scattering matrix is projectile, second index is target */
                        if (result != 0) {
                            printf ("Error: insufficient memory for target-target scattering matrices! %i \n",
                                    result);
                            return result;
                        }
                    }
                }
            }
        }

        /* For each element, possibly create arrays for storing leaving recoils */
        if (store_exiting_recoils == 1) {
            /* create arrays first: */
            for (i=0; i<number_of_materials; i++) {
                list_of_materials[i].elemental_leaving_recoils =
                    malloc (sizeof (struct transmitted_ion*) * list_of_materials[i].element_count);
                list_of_materials[i].leaving_recoils_pointer =
                    malloc (sizeof (int) * list_of_materials[i].element_count);
                for (j=0; j<list_of_materials[i].element_count; j++) {  /* go through elements */
                    (list_of_materials[i].elemental_leaving_recoils)[j] =
                        malloc(store_exiting_limit * sizeof (struct transmitted_ion));
                    if ((list_of_materials[i].elemental_leaving_recoils)[j] == NULL) {
                        printf ("Error: Insufficient memory for arrays of leaving recoils!\n");
                        return -2001;
                    }
                    (list_of_materials[i].leaving_recoils_pointer)[j] = 0;
                }
            }
        }
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=0) printf (" finished\n"); fflush(stdout);
        /* MPI=============================================== */
#else
        if (print_level >= 0) printf (" finished\n"); fflush(stdout);
#endif
    }
    else {  /* only estimate memory usage */
#ifndef MPI_PRALLEL
        /* MPI=============================================== */
        if (store_exiting_recoils == 1) {
            lui_temp = 0;
            for (i=0; i<number_of_materials; i++) {
                lui_temp += (sizeof (struct transmitted_ion*) + sizeof (int) + store_exiting_limit
                          *  sizeof (struct transmitted_ion)) * list_of_materials[i].element_count;
            }
            mem_usage += lui_temp;
            if (mem_usage_details == 1)
                printf ("MEMORY Leaving recoils arrays:      %lu bytes\n", lui_temp);
        }
        number_exist_element = count_existing_elements (existing_elements);  /* (1) - (7) - (5) */
        lui_temp = number_exist_element * (number_of_materials - 1) * MAX_STOPPING_ENTRIES * sizeof (float);
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Stopping tables:             %lu bytes\n", lui_temp);
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Straggling tables:           %lu bytes\n", lui_temp);
        lui_temp += number_exist_element * 2 * (DIME * DIMS * sizeof (float));
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Ion scattering matrices:     %li bytes\n", lui_temp);
        lui_temp += number_exist_element * number_exist_element * 2 * (DIME * DIMS * sizeof (float));
        mem_usage += lui_temp;
        if (mem_usage_details == 1)
            printf ("MEMORY Recoil scattering matrices:  %li bytes\n", lui_temp);
        /* MPI=============================================== */
#endif
    }

    return 0;
}

/*=============================================================================
  Function Name : read_materials_data_block
  Description   : Needs to be called from the init file reader while the materials
                  input file is read and the definition of a new material begins.

  Inputs  : char* material_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int read_materials_data_block (char *material_name) {

    /* check, if there is still space left in the material list */
    if (number_of_materials >= MAX_NO_MATERIALS) {  /* material index is full */
        printf ("Warning! Maximum number of materials (%i) exceeded!\n", MAX_NO_MATERIALS);
        return -2002;
    }

    /* make sure, that material-name is not to long, i.e. cut it if it is */
    if (strlen (material_name) >= 25) material_name[24] = '\0';

    strcpy (list_of_materials[number_of_materials].name, material_name);  /* store Name */
    list_of_materials[number_of_materials].element_count   = 0;
    list_of_materials[number_of_materials].density         = 0;
    list_of_materials[number_of_materials].ion_surf_energy = -0.001;

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=1) printf ("\nNew material declared:\t%s\n", material_name);
    /* MPI=============================================== */
#else
    if (print_level >= 1) printf ("\nNew material declared:\t%s\n", material_name);
#endif

    number_of_materials++;  /* increase index */

    return 0;
}

/*=============================================================================
  Function Name : read_materials_data
  Description   : Needs to be called from the init file reader while the materials
                  input file is read.

  Inputs  :
            char* par_name
            char* par_value
  Outputs : no.

  Notes : no.
=============================================================================*/
int read_materials_data (char *par_name, char *par_value) {

    /* check some limits */
    if (number_of_materials >= MAX_NO_MATERIALS) return -2003;  /* too large */

    /* Compare the parameter name to known parameters and then read in corresponding value */

    if (strcmp (par_name, "element_count") == 0) {  /* read number of element into the
                                                     current material */
        sscanf (par_value,"%i", &(list_of_materials[number_of_materials-1].element_count));
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=1)
            printf ("Number of elements:\t%i\n",
                    list_of_materials[number_of_materials-1].element_count);
        /* MPI=============================================== */
#else
        if (print_level >= 1)
            printf ("Number of elements:\t%i\n",
                    list_of_materials[number_of_materials-1].element_count);
#endif
    }
    if (strcmp (par_name, "density") == 0) {  /*  */
        sscanf (par_value, "%f", &(list_of_materials[number_of_materials-1].density));
        if (isnan (list_of_materials[number_of_materials-1].density))
            list_of_materials[number_of_materials-1].density = 0.0;
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=1)
            printf ("Density:\t\t%g\n", list_of_materials[number_of_materials-1].density);
        /* MPI=============================================== */
#else
        if (print_level >= 1)
            printf ("Density:\t\t%g\n", list_of_materials[number_of_materials-1].density);
#endif
    }
    if (strcmp (par_name,"elements_Z") == 0)  /* read list of integer Z values*/
        if (make_int_array (par_value, MAX_EL_PER_MAT,
            list_of_materials[number_of_materials-1].elements_Z) != 0) return -2004;
    if (strcmp (par_name, "elements_M") == 0)  /* read list of mass values*/
        if (make_float_array (par_value, MAX_EL_PER_MAT,
            list_of_materials[number_of_materials-1].elements_M) != 0) return -2005;
    if (strcmp (par_name, "elements_conc") == 0)  /* read list of concentrations */
        if (make_float_array(par_value,MAX_EL_PER_MAT,
            list_of_materials[number_of_materials-1].elements_conc) != 0) return -2006;
    if (strcmp (par_name, "elements_disp_energy") == 0)  /* read list of displacement energies */
        if (make_float_array (par_value, MAX_EL_PER_MAT,
            list_of_materials[number_of_materials-1].elements_disp_energy) != 0) return -2007;
    if (strcmp (par_name, "elements_latt_energy") == 0)  /* read list of lattice energies */
        if (make_float_array (par_value, MAX_EL_PER_MAT,
            list_of_materials[number_of_materials-1].elements_latt_energy) != 0) return -2008;
    if (strcmp (par_name, "elements_surf_energy") == 0)  /* read list of surface binding energies */
        if (make_float_array(par_value, MAX_EL_PER_MAT,
            list_of_materials[number_of_materials-1].elements_surf_energy) != 0) return -2009;
    if (strcmp (par_name, "ion_surf_energy") == 0) {  /* read ion surface binding energy */
        sscanf (par_value,"%f", &(list_of_materials[number_of_materials-1].ion_surf_energy));
#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node==ROOT && print_level>=1)
            printf ("Ion SBE:\t\t%g\n", list_of_materials[number_of_materials-1].ion_surf_energy);
        /* MPI=============================================== */
#else
        if (print_level >= 1)
            printf ("Ion SBE:\t\t%g\n", list_of_materials[number_of_materials-1].ion_surf_energy);
#endif
    }

    return 0;
}

/*=============================================================================
  Function Name : prepare_stopping_tables
  Description   : Read stopping data from file and fill arrays etc.
                  Material version.

  Inputs  : no.
  Outputs : no.

  Notes :
      Function call:
        fileio - int read_float_block (char* file_name, int offset, int count,
                                       float* array);
=============================================================================*/
int prepare_stopping_tables (void) {
    int i, j, k, l, result;
    char StoppingFileName[1000];
    float flt_temp[DIMD+1];

    /* go through materials and create stopping tables for all existing elements */
    for (i=0; i<number_of_materials; i++) {
        /* pointers to stopping arrays */
        list_of_materials[i].stopping_ZE = (float**) malloc (sizeof (float*) * MAX_ELEMENT_NO);
        if (list_of_materials[i].stopping_ZE == NULL) return -2010;

        /* go through possible elements which might be slowed down in this materials */
        for (j=0; j<MAX_ELEMENT_NO; j++) {  /* goes through all possible projectiles */
            if (existing_elements[j] == 1) {  /* ok, element might occur, calculate */
                /* allocate memory for stopping table */
                list_of_materials[i].stopping_ZE[j] = (float*) calloc (MAX_STOPPING_ENTRIES, sizeof (float));
                if (list_of_materials[i].stopping_ZE[j] == NULL) return -2011;

                /* OK, now we can load the stopping table from the file */
                /* the following code to do this is adapted from the corteo code */
                sprintf (StoppingFileName, "data/%u.asp", j);  /* filename where stopping data are tabulated */

                /* for all elements occurring in the current material, we need to load
                   the stopping data and then apply a rule for stopping in compounds */
                /* go through elements of target material */
                for (k=0; k<list_of_materials[i].element_count; k++) {
                    result = read_float_block (StoppingFileName,
                                               (DIMD + 1) * (list_of_materials[i].elements_Z[k] - 1),
                                                DIMD + 1, flt_temp);
                    if (result != 0) return -2012;  /* could not read float data block from file */

                    /* each record must end with the element number as control data */
                    if (flt_temp[DIMD] != list_of_materials[i].elements_Z[k]) return -2013;
                    for (l=0; l<DIMD; l++) {
                        /* add stopping for current target element to stopping of current target material */
                        /* electronic stopping of index k; factor 10 is conversion from eV/(1E15 at/cm^2)
                           to eV/(at/A^2), density is in at/A^3 (=10^24 at/cm^3) */
                        (list_of_materials[i].stopping_ZE[j])[l] += flt_temp[l]
                                        * list_of_materials[i].elements_conc[k] * 10.0f
                                        * list_of_materials[i].density * 1e-24;  /* must be in at/A^3 CHECK */
                                        /* * compoundCorr; */ /* Bragg's rule */

                        /* TODO: the compound correction needs to be added here !!! */
                        /* following copied from corteo:
                           compound correction according to Zeigler & Manoyan NIMB35(1988)215, Eq.(16)
                           (error in Eq. 14)
                           if(compoundCorr!=1.0f)
                           for(k=0; k<DIMD; k++) {
                               f = d2f (1.0 / (1.0 + exp (1.48 * (sqrt (2. * Dval(k) / projectileMass / 25e3) - 7.0))));
                               spp[k] *= f * (compoundCorr - 1.0f) + 1.0f;
                           }
                        */
                        /* unit here is eV/nm */
                    }  /* end of loop through index stopping values */
                }  /* end of loop through elements of target material */
            }
            else {  /* if element does not occur, create NULL pointer */
                list_of_materials[i].stopping_ZE[j] = NULL;
            }
        }  /* end of loop through possible projectile elements */
    }  /* end of loop through target materials */

    return 0;
}

/*=============================================================================
  Function Name : prepare_straggling_tables
  Description   : Create straggling tables etc.
                  Material version!
                  model:
                      0: no straggling
                      1: Bohr straggling
                      2: Chu correction        PRA  13 (1976) 2057
                      3: Chu + Yang correction NIMB 61 (1991) 149
                  The function load_Chu_straggling_values() must have been
                  called before.

  Inputs  : int model
  Outputs : no.

  Notes :
      Function call:
        fromcorteo - float Dval (unsigned long index);
                     unsigned long D_index (float val);
                     loat d2f (double val);
                     loat sqrtdf (double val).
=============================================================================*/
int prepare_straggling_tables (int model) {
    int i, Z, k, l;  /* projectile's Z */
    unsigned long ii;

    double straggling;
    double stragg_element;
    double stopping;
    double energy;
    double MEV_energy_amu;
    double chargestate2;  /* the effective charge state of current projectile */
    double mass;  /* mass of most abundant isotope of projectile Z */
    int    target_Z;
    double omega_Bohr2;
    double Chu_factor;  /* Chu's correction factor for the Bohr straggling */
    double Yang, epsilon, gamma;
    double C1, C2, C3, C4, B1, B2, B3, B4;

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==0 && print_level>1) printf ("Straggling model %i\n", model);
    /* MPI=============================================== */
#else
    if (print_level > 1) printf ("Straggling model %i\n", model);
#endif

    /* go through materials and create straggling tables for all existing elements */
    for (i=0; i<number_of_materials; i++) {
        /* pointers to straggling arrays */
        list_of_materials[i].straggling_ZE = (float**) malloc (sizeof(float*) * MAX_ELEMENT_NO);
        if (list_of_materials[i].straggling_ZE == NULL) return -2014; /* cannot allocate memory */

        for (Z=0; Z<MAX_ELEMENT_NO; Z++) {  /* loop through all possible projectiles */
            if (existing_elements[Z] == 1) {  /* ok, element might occur as projectile, calculate */
                /* allocate memory for straggling table of projectile j in material i */
                list_of_materials[i].straggling_ZE[Z] = (float*) calloc (MAX_STOPPING_ENTRIES, sizeof (float));
                if (list_of_materials[i].straggling_ZE[Z] == NULL) return -2015;

                mass = most_abundant_isotope[Z];

                /* calculation of straggling (similar to corteo code): */
                for(k=0; k<DIMD; k++) {  /* go through table that has to be filled */
                    straggling = 0;

                    for (l=0; l<list_of_materials[i].element_count; l++){  /* All elements of current target material */
                        target_Z = list_of_materials[i].elements_Z[l];
                        stopping = (list_of_materials[i].stopping_ZE[Z])[k];
                        energy = D_val (k);  /* energy corresponding to index k */

                        /* the effective charge state is obtained from comparing stopping of the ion
                           and the hydrogen:
                        chargestateqaured = stopping(H)/stopping(ion) Z_ion^2 for stopping at same speed */
                        ii = D_index (d2f (energy / mass));  /* index of the velocity which is
                                                               proton energy of same velocity */
                        chargestate2 = stopping / (((list_of_materials[i].stopping_ZE[1])[ii]) * Z * Z);

                        /* start by calculating squared Bohr straggling (all other models need this anyway) */
                        omega_Bohr2 = 4.0 * PI * Z * Z * target_Z *  E2 * E2
                                   * list_of_materials[i].density * 1e-24  /* must be in at/A^3 CHECK */;

                        /* calculate the Chu correction factor. the formula needs energy[MeV]/mass: */
                        MEV_energy_amu = energy * 1e-6 / mass;
                        Chu_factor = 1.0 / ( 1.0 + chu_values[Z][0] * pow (MEV_energy_amu, chu_values[Z][1])
                                   + chu_values[Z][2] * pow (MEV_energy_amu, chu_values[Z][3]));

                        /* to calculate Yang's extra straggling contribution caused by charge state
                           fluctuations, we need his Gamma and his epsilon (eq.6-8 from the paper): */
                        /* For hydrogen we need the B and for other projectile the C constants: */
                        if (Z == 1) {
                            B1 = 0.1955; B2 = 0.6941; B3 = 2.522; B4 = 1.040;
                            gamma = B3 * (1.0- exp (-B4 * MEV_energy_amu));
                            Yang  = B1 * gamma / (pow ((MEV_energy_amu-B2), 2.0) + gamma * gamma);
                        }
                        else {
                            C1 = 1.273e-2; C2 = 3.458e-2; C3 = 0.3931; C4 = 3.812;  /* solid targets */
                            epsilon = MEV_energy_amu *  pow (Z, -1.5)  * pow (target_Z, -0.5);  /* solid targets */
                            gamma = C3 * (1.0 - exp (-C4 * epsilon));
                            Yang  = (pow (Z, 1.333333333333) / pow (target_Z, 0.33333333333)) * C1 * gamma
                                  / (pow ((epsilon - C2), 2.0) + gamma * gamma);
                        }

                        /* Now we have all ingredients for any of the straggling models. We could have saved some
                           calculations by checking the model first, but well... we'll probably use Yang's model
                           in most cases */
                        switch (model) {
                        case 0 :  /* no straggling */
                            stragg_element = 0.;
                            break;
                        case 1 :  /* Bohr */
                            stragg_element = omega_Bohr2;
                            break;
                        case 2 :  /* Chu */
                            stragg_element = omega_Bohr2 * Chu_factor;
                            break;
                        case 3 :  /* Chu + Yang correction */
                            stragg_element = omega_Bohr2 * (chargestate2 * Chu_factor + Yang);
                            break;
                        default :
                            stragg_element = 0;
                        }

                        /* Now we know the straggling for each target element in the material
                           we can add them up using Bragg's rule of additivity */
                        straggling += stragg_element * list_of_materials[i].elements_conc[l];
                    }  /* end of loop through elements in current target material, l */

                    /* store the straggling in its table: */
                    (list_of_materials[i].straggling_ZE[Z])[k] = sqrtdf (straggling) * sqrtdf (2.0);
                }  /* end of loop through entries in straggling table, k */
            }  /* end the check if element might occur as projectile */
        }  /* end of possible projectiles loop */
    }  /* end of target material loop */

    return 0;
}

/*=============================================================================
  Function Name : prepare_KP_tables
  Description   : Modified Kinchin-Pease model, Material version.

  Inputs  : no.
  Outputs : no.

  Notes :
      Added on Aug. 05, 2014.
=============================================================================*/
int prepare_KP_tables1 (void) {
    int i, j;
    double target_Z, target_M, target_Ed, energy, k_d, g_ed, e_d, E_v, E_c;

    for (i=0; i<number_of_materials; i++) {
        //list_of_materials[i].Nd = (float*) malloc (sizeof (float) * MAX_STOPPING_ENTRIES);
        list_of_materials[i].Nd = (float*) calloc (MAX_STOPPING_ENTRIES, sizeof (float));
        if (list_of_materials[i].Nd == NULL) return -2016;

        target_Z = list_of_materials[i].mean_Z;
        target_M = list_of_materials[i].mean_M;
        target_Ed = list_of_materials[i].mean_Ed;
        E_c = 2.5 * target_Ed;

        for(j=0; j<DIMD; j++) {  /* go through table that has to be filled */
            energy = D_val (j);  /* energy corresponding to index j */

            k_d = 0.133745 * pow (target_Z, 2.0 / 3.0) / pow (target_M, 0.5);

            e_d = 0.0115 * pow (target_Z, -7.0 / 3.0) * energy;
            g_ed = 3.4008 * pow (e_d, 1.0 / 6.0) + 0.40244 * pow (e_d, 0.75) + e_d;

            E_v = energy / (1.0 + k_d * g_ed);

            if (E_v < target_Ed) {
                list_of_materials[i].Nd[j] = 0.0;
            }
            else if (E_v >= target_Ed && E_v < E_c) {
                list_of_materials[i].Nd[j] = 1.0;
            }
            else if (E_v >= E_c) {
                list_of_materials[i].Nd[j] = E_v / E_c;
            }
            //printf ("0000: %f\n", list_of_materials[i].Nd[j]);
        }
    }

    return 0;
}

/*=============================================================================
  Function Name : prepare_KP_tables
  Description   : Modified Kinchin-Pease model, Material version.

  Inputs  : no.
  Outputs : no.

  Notes :
      Wrong!!!
      Added on Aug. 05, 2014.
=============================================================================*/
//int prepare_KP_tables2 (void) {
//    int i, j, k;
//    double target_Z, target_M, target_Ed, energy, k_d, g_ed, e_d, E_v, E_c;

//    for (i=0; i<number_of_materials; i++) {
//        list_of_materials[i].Nd = (float*) malloc (sizeof (float) * MAX_STOPPING_ENTRIES);
//        if (list_of_materials[i].Nd == NULL) return -2111;
//
//        for (j=0; j<list_of_materials[i].element_count; j++) {
//            target_Z = list_of_materials[i].elements_Z[j];
//            target_M = list_of_materials[i].elements_M[j];
//            target_Ed = list_of_materials[i].elements_disp_energy[j];
//            E_c = 2.5 * target_Ed;

//            for(k=0; k<DIMD; k++) {  /* go through table that has to be filled */
//                energy = D_val (k);  /* energy corresponding to index j */

//                k_d = 0.133745 * pow (target_Z, 2.0 / 3.0) / pow (target_M, 0.5);

//                e_d = 0.0115 * pow (target_Z, -7.0 / 3.0) * energy;
//                g_ed = 3.4008 * pow (e_d, 1.0 / 6.0) + 0.40244 * pow (e_d, 0.75) + e_d;

//                E_v = energy / (1.0 + k_d * g_ed);

//                if (E_v < target_Ed) {
//                    list_of_materials[i].Nd[k] += 0.0;
//                }
//                else if (E_v >= target_Ed && E_v < E_c) {
//                    list_of_materials[i].Nd[k] += list_of_materials[i].elements_conc[j];
//                }
//                else if (E_v >= E_c) {
//                    list_of_materials[i].Nd[k] += (E_v / E_c) * list_of_materials[i].elements_conc[j];
//                }
                //printf ("0000: %f\n", list_of_materials[i].Nd[k]);
//            }
//        }
//    }

//    return 0;
//}

/*=============================================================================
  Function Name : prepare_KP_tables
  Description   : Modified Kinchin-Pease model, Material version.

  Inputs  : no.
  Outputs : no.

  Notes :
      Z_1 and Z_2, M_1 and M_2.
      Added on Aug. 13, 2014.
=============================================================================*/
int prepare_KP_tables2 (void) {
    int i, j, k, l;
    double project_Z, project_M, target_Z, target_M, target_Ed, energy, k_d,
           g_ed, aa, e_d, E_v, E_c;

    for (i=0; i<number_of_materials; i++) {
        list_of_materials[i].Nd_Z = (float**) calloc (list_of_materials[i].element_count, sizeof (float*));
        //list_of_materials[i].Nd_Z = (float**) malloc (sizeof (float*) * list_of_materials[i].element_count);
        if (list_of_materials[i].Nd_Z == NULL) return -2017;

        for (j=0; j<list_of_materials[i].element_count; j++) {
            project_Z = list_of_materials[i].elements_Z[j];
            project_M = list_of_materials[i].elements_M[j];

            list_of_materials[i].Nd_Z[j] = (float*) malloc (sizeof (float) * MAX_STOPPING_ENTRIES);
            if (list_of_materials[i].Nd_Z == NULL) return -2018;


            for (k=0; k<list_of_materials[i].element_count; k++) {
                target_Z = list_of_materials[i].elements_Z[k];
                target_M = list_of_materials[i].elements_M[k];
                target_Ed = list_of_materials[i].elements_disp_energy[k];
                E_c = 2.5 * target_Ed;

                for(l=0; l<DIMD; l++) {  /* go through table that has to be filled */
                    energy = D_val (l);  /* energy corresponding to index j */

                    k_d = 0.133745 * pow (project_Z, 2.0 / 3.0) / pow (project_M, 0.5);

                    aa = 0.0325 * pow (pow (project_Z, 2.0 / 3.0) + pow (target_Z, 2.0 / 3.0), -0.5);
                    e_d = target_M * energy / (project_M + target_M) * aa / (project_Z * target_Z);
                    g_ed = 3.4008 * pow (e_d, 1.0 / 6.0) + 0.40244 * pow (e_d, 0.75) + e_d;

                    E_v = energy / (1.0 + k_d * g_ed);

                    if (E_v < target_Ed) {
                        list_of_materials[i].Nd_Z[j][l] += 0.0;
                    }
                    else if (E_v >= target_Ed && E_v < E_c) {
                        list_of_materials[i].Nd_Z[j][l] += list_of_materials[i].elements_conc[k];
                    }
                    else if (E_v >= E_c) {
                        list_of_materials[i].Nd_Z[j][l] += (E_v / E_c) * list_of_materials[i].elements_conc[k];
                    }
                    //printf ("0000: %f\n", list_of_materials[i].Nd_Z[j][l]);
                }
            }
        }
    }

    return 0;
}

/*=============================================================================
  Function Name : count_existing_elements
  Description   : Returns the number of ones in the provided array.

  Inputs  : int* elementarray
  Outputs : no.

  Notes : no.
=============================================================================*/
int count_existing_elements (int *element_array) {
    int result = 0;
    int i;

    for (i=0; i<MAX_ELEMENT_NO; i++) result += element_array[i];

    return result;
}

/*=============================================================================
  Function Name : convert_material_to_element
  Description   : Reads in the "standard" config file for material-based target
                  definition and creates files for element based target definition.

  Inputs  : char* output_file
  Outputs : no.

  Notes :
      This function only works in the non-dynamic compilation, because this way
      it is easier to write, since the standard read-procedures can be used.

      Function call:
        target - int get_target_XYZ (int index, int* x, int* y, int* z);
        fileio - int combine_files (int count, ...).
=============================================================================*/
int convert_material_to_element (char *output_file) {
    int i, j, k;
    int x, y, z;                /* integer cell indices in each direction */
    struct material* cell_mat;  /* pointer to material of current cell */
    FILE*  elem_fp;             /* points to element output file */
    FILE*  comp_fp;             /* points to new composition file */
    char   str_temp[MAX_FILENAME_LENGTH];  /* for temporary filenames */
    int    result;

    float mean_disp;  /* mean energies are calculated and used for the ion...
                         its simple and may not be correct, */
    float mean_latt;  /* but we need some values. Can be changed later by the user */
    float mean_surf;
    float mean_mass;
    int   elem_count, elem_count2;

    float* compVector;      /* for writing the composition file. Has as many entries
                               as different elements occur */
    int* mat_elem_vector_p;  /* for each material and element, we will store here the
                                index of the composition vector, which corresponds to
                                the Z of the element. We need this for constructing the
                                composition file. It is not part of the material structure,
                                because it is not used for normal operation of im3d */
    char str_temp2[1024];

    mean_disp   = .0f;
    mean_latt   = .0f;
    mean_surf   = .0f;
    mean_mass   = .0f;
    elem_count  = 0;
    elem_count2 = 0;

    printf ("\nAttempting to convert material based target definition to element based target composition... \n\n");

    mat_elem_vector_p = (int*) malloc (MAX_NO_MATERIALS * MAX_EL_PER_MAT * sizeof (int));
    if (mat_elem_vector_p == NULL) {
        printf ("Error: insufficient memory!\n");
        return -2019;
    }

    if (single_input_file == 1) {  /* OK, create also a single output file. So we need a temporary
                      element file to write to, and a temporary composition file */
        strcpy (ElementsFileName, "temp_elementfile.im3d");
        strcpy (TargetCompositionFileName, "temp_compfile.im3d");
        strcpy (ConversionFileName, output_file);
    }
    else {  /* use given output file name to put element definition to */
        strcpy (ElementsFileName, output_file);
        sprintf (str_temp, "%s.e", TargetCompositionFileName);  /* and append .e to new composition file */
        strcpy (TargetCompositionFileName, str_temp);
    }
    elem_fp = fopen (ElementsFileName, "w");
    if (elem_fp == NULL) {
        printf ("Error. Cannot open file %s for writing.\n", ElementsFileName);
        return -2020;
    }

    /* create element file */
    if (conv_create_separate_elements == 1) {  /* create separate elements */
        if (print_level > -2)
            printf ("Parsing elements (create separate elements for each material) ...\n");
        fprintf (elem_fp, "ElementCount=%i\n\n", number_of_all_target_elements);
        for (i=0; i<number_of_materials; i++) {
            for (j=0; j<list_of_materials[i].element_count; j++) {
                if (print_level > 0)
                    printf ("Element %s found in material %s. Storing data.\n",
                        atomic_names[list_of_materials[i].elements_Z[j]], list_of_materials[i].name);
                fprintf (elem_fp, "[%s_from_%s]\n",
                         atomic_names[list_of_materials[i].elements_Z[j]], list_of_materials[i].name);
                fprintf (elem_fp, "Z=%i\n", list_of_materials[i].elements_Z[j]);
                fprintf (elem_fp, "M=%g\n", list_of_materials[i].elements_M[j]);
                fprintf (elem_fp, "DispEnergy=%g\n", list_of_materials[i].elements_disp_energy[j]);
                fprintf (elem_fp, "LattEnergy=%g\n", list_of_materials[i].elements_latt_energy[j]);
                fprintf (elem_fp, "SurfEnergy=%g\n\n", list_of_materials[i].elements_surf_energy[j]);
                mean_disp += list_of_materials[i].elements_disp_energy[j];  /* add stuff up to calc mean values */
                mean_latt += list_of_materials[i].elements_latt_energy[j];
                mean_surf += list_of_materials[i].elements_surf_energy[j];
                mat_elem_vector_p[i*MAX_EL_PER_MAT+j] = elem_count+1;  /* note, which element number this will get later */
                if (print_level > 2)
                    printf ("New element number %i assigned to mat. %i, elem. %i\n", elem_count + 1, i, j);
                elem_count ++;
            }
        }
        /* calculate mean values */
        mean_disp /= (float) elem_count;
        mean_latt /= (float) elem_count;
        mean_surf /= (float) elem_count;
        if (print_level > -2) printf ("Finished parsing elements.\n");
    }
    else {  /* each element to appear only once */
        if (print_level > -2) printf ("Parsing elements (list each elements not more than once) ...\n");
        /* first: count element to put into file: */
        for (i=1; i<MAX_ELEMENT_NO; i++) {
            if (existing_elements[i] == 1) {  /* elements exists! */
                if (((i==1)     && (hydrogen_in_target==1)) ||
                    ((i>1)      && (i != ion_Z)) ||
                    ((i==ion_Z) && (ionZ_in_target==1))) {  /* OK, element really exists in target,
                                                    not just additionally added ion or hydrogen */
                    elem_count ++;
                }
            }
        }
        fprintf (elem_fp, "ElementCount=%i\n\n", elem_count);
        elem_count = 0;  /* restart counting !*/
        /* then work on elements */
        for (i=1; i<MAX_ELEMENT_NO; i++) {
            if (existing_elements[i] == 1) {  /* elements exists! */
                if (((i==1)     && (hydrogen_in_target==1)) ||
                    ((i>1)      && (i!=ion_Z)) ||
                    ((i==ion_Z) && (ionZ_in_target==1))) {  /* OK, element really exists in target,
                                                    not just additionally added ion or hydrogen */
                    if (print_level > 0) printf ("Element %s found. Storing data.\n", atomic_names[i]);
                    /* Now search in which materials the element appears and obtain mean values
                       for mass, energies... : */
                    mean_disp = 0; mean_latt = 0; mean_surf = 0; mean_mass = 0; elem_count2 = 0;
                    for (k=0; k<number_of_materials; k++) {
                        for (j=0; j<list_of_materials[k].element_count; j++) {
                            if (list_of_materials[k].elements_Z[j] == i) {  /* element found in material */
                                /* note, which element number this will get later */
                                mat_elem_vector_p[k*MAX_EL_PER_MAT+j] = elem_count + 1;
                                if (print_level>2)
                                    printf ("New element number %i assigned to mat. %i, elem. %i. Absolute element: %i\n",
                                        elem_count + 1, k, j, i);
                                /* add stuff up to calc mean values */
                                mean_disp += list_of_materials[k].elements_disp_energy[j];
                                mean_latt += list_of_materials[k].elements_latt_energy[j];
                                mean_surf += list_of_materials[k].elements_surf_energy[j];
                                mean_mass += list_of_materials[k].elements_M[j];
                                elem_count2++;
                            }
                        }
                    }
                    mean_disp /= (float) elem_count2;
                    mean_latt /= (float) elem_count2;
                    mean_surf /= (float) elem_count2;
                    mean_mass /= (float) elem_count2;
                    /* Ok, write element info: */
                    fprintf (elem_fp, "[%s]\n", atomic_names[i]);
                    fprintf (elem_fp, "Z=%i\n", i);
                    fprintf (elem_fp, "M=%g\n", mean_mass);
                    fprintf (elem_fp, "DispEnergy=%g\n", mean_disp);
                    fprintf (elem_fp, "LattEnergy=%g\n", mean_latt);
                    fprintf (elem_fp, "SurfEnergy=%g\n\n", mean_surf);
                    elem_count ++;
                }
                else {  /* element is in list, but does not really exist in target, ignore */
                }
            }
        }
        if (print_level > -2) printf ("Finished parsing elements.\n");
           /* For now, the ion's displacement, lattice and surface energy are simply
              determined by the last element that appeared.
              The values should anyway better be set by the user later! */
    }

    fprintf (elem_fp, "[ion]\n");
    fprintf (elem_fp, "DispEnergy=%g\n", mean_disp);
    fprintf (elem_fp, "LattEnergy=%g\n", mean_latt);
    fprintf (elem_fp, "SurfEnergy=%g\n", mean_surf);

    fclose (elem_fp);

    /* OK, element file has been written, new create composition file */
    if (print_level > -2) printf ("New element file has been created.\n");

    compVector = (float*) malloc ((elem_count + 1) * sizeof (float));
    if (compVector == NULL) {
        printf ("Error: insufficient memory!\n");
        return -2021;
    }

    comp_fp = fopen (TargetCompositionFileName, "w");
    if(comp_fp==NULL) {
        printf ("Error. Cannot open file %s for writing.\n", TargetCompositionFileName);
        return -2022;
    }

    if (print_level > -2)
        printf ("Creating new composition file %s ...\n", TargetCompositionFileName);
    for (i=0; i<cell_count; i++) {  /* cycle through all cells and write entries to file */
        get_target_XYZ (i, &x, &y, &z);  /* get coords of cell */
        cell_mat = &(list_of_materials[target_composition[i]]);  /* get material of cell */
        fprintf (comp_fp, "%i\t%i\t%i\t", x, y, z);  /* print coords to file */
        fprintf (comp_fp, "%g", cell_mat->density);  /* print density to file */

        /* now construct composition vector for current cell: */
        for (j=0; j<=elem_count; j++) compVector[j] = 0;  /* init with zeros */
        for (j=0; j<cell_mat->element_count; j++)  /* cycle through elements of this material and get the concentration */
            compVector[mat_elem_vector_p[target_composition[i]*MAX_EL_PER_MAT+j]] = cell_mat->elements_conc[j];

        /* now, write the composition vector into the composition file */
        for (j=0; j<=elem_count; j++) fprintf (comp_fp, "\t%g", compVector[j]);
        fprintf (comp_fp, "\n");
    }
    fclose (comp_fp);
    if (print_level > -2)
        printf ("New composition file %s has been created.\n", TargetCompositionFileName);

    if (single_input_file == 1) {  /* OK, create a single output file */
        if (print_level > -2)
            printf ("Combining input files to create single project file %s... \n", ConversionFileName);
        sprintf (str_temp2, "ElementsFileName=%s\n\n#<<<BEGIN STRUCTUREFILE", ElementsFileName);
        result = combine_files (9, ConversionFileName, "#<<<BEGIN CONFIGFILE", "temp_configfile.im3d",
                 str_temp2, "temp_structfile.im3d", "#<<<BEGIN ELEMFILE", ElementsFileName,
                 "#<<<BEGIN COMPFILE", TargetCompositionFileName);
        if (result != 0) {
            printf ("Error %i. Cannot create combined output file.\n", result);
            return result;
        }
        if (print_level > -2)
            printf ("Combined input files %s has been created.\n", ConversionFileName);
    }

    if (print_level > -2)
        printf("\nMaterial based input file successfully converted to element based input file. \n\n");

    return 0;
}
