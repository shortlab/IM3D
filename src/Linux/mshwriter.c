/******************************************************************************
  Module Name : mshwriter.h
  Module Date : 05/02/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Output files with MSH ASCII/binary file format.

  Others :
      Error numbers in this module 6000-6999.
      Gmsh software's native "MSH" file format.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "mshwriter.h"

/*=============================================================================
  Function Name : store_results_msh
  Description   : Store the results of the simulation (arrays with distribution
                  of implanted ions, defects etc.) in .msh format.
                  This is the non-dynamic version of the function for
                  material-based output.

  Inputs  :
            char* base_name
  Outputs : no.

  Notes : ygli
=============================================================================*/
int store_results_msh (char *base_name) {
    int  base_name_length;
    char *str_temp;
    char *str_temp2;
    int  i, j;
    FILE *fp;

    base_name_length = strlen (base_name);
    str_temp  = (char*) malloc (sizeof (char) * (base_name_length + 100));
    str_temp2 = (char*) malloc (sizeof (char) * (255));
    strncpy (str_temp, base_name, base_name_length + 1);

    /* Arrays to be stored:
       - the array with concentration of implanted ions
       - sum of displacements, replacement and so on in each cell
       - in case the program is ever made dynamic: the final concentration
       - for each element of each material: vacancies and interstitial arrays
       - transmitted ions
       - sputtered atoms in each direction
    */

    /* concentration of implanted ions: */
    strcat (str_temp, "ions.total.msh");
    if (print_level >= 2) printf ("Storing implanted ions to:      %s\n", str_temp);
    write_int_array_to_msh_file (str_temp, 0, target_implanted_ions, 0, cell_count, 0);
    str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */

    if (do_not_store_damage == 0) {
        /* Part of ions that replaced identical target atoms */
        strcat (str_temp, "ions.replacements.msh");
        if (print_level >= 2) printf ("Storing implanted replacing ions to:   %s\n", str_temp);
        write_int_array_to_msh_file (str_temp, 0, target_replacing_ions, 0, cell_count, 0);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    }

    sum_up_material_arrays ();  /* sum up the individual element-dependent arrays */

    /* deposited energy */
    if (store_energy_deposit == 1) {
        strcat (str_temp, "energy.deposit.msh");
        if (print_level >= 2) printf ("Storing energy deposited to electrons and phonons:  %s\n", str_temp);
        write_double_array_to_msh_file (str_temp, target_energy_electrons, target_energy_phonons, cell_count);
        str_temp[base_name_length] = '\0';
    }

    /* store single ion sputter yields */
    if (single_ion_sputter_yields == 1) {
        strcat (str_temp, "single_ion_sputter_yields");
        if (print_level >= 2) printf ("Storing single ion sputter yields to:   %s\n", str_temp);
        fp = fopen (str_temp, "w");
        if (fp == NULL) {
            printf("Error: cannot store histogram.\n");
            return -6000;
        }
        for (i=0; i<=MAX_SPUTTERED; i++) fprintf (fp, "%i\t%i\n", i, sputter_yield_histogram[i]);
        fclose (fp);
        str_temp[base_name_length] = '\0';
    }

    /* sum of atoms leaving the target */
    if (detailed_sputtering == 1) {
        /* particles leaving the sample (sputtered or implanted deeper): */
        strcat (str_temp, "leaving_directions.sum");
        if (print_level >= 2) printf (" Storing sum of leaving atoms to:       %s\n", str_temp);
        sprintf (str_temp2, "(+x,+y,+z)\t%i\n(+x,-y,+z)\t%i\n(+x,-y,-z)\t%i\n(+x,+y,-z)\t%i\n(-x,+y,+z)\t%i\n(-x,-y,+z)\t%i\n(-x,-y,-z)\t%i\n(-x,+y,-z)\t%i\n",
                 total_sputter_counter[0], total_sputter_counter[1], total_sputter_counter[2], total_sputter_counter[3],
                 total_sputter_counter[4], total_sputter_counter[5], total_sputter_counter[6], total_sputter_counter[7]);
        write_string_to_file (str_temp, str_temp2);
        str_temp [base_name_length] = '\0';

        /* ions leaving the target */
        strcat (str_temp, "leaving_directions.ions");
        if (print_level >= 2) printf (" Storing sum of leaving ions to:        %s\n", str_temp);
        sprintf (str_temp2, "(+x,+y,+z)\t%i\n(+x,-y,+z)\t%i\n(+x,-y,-z)\t%i\n(+x,+y,-z)\t%i\n(-x,+y,+z)\t%i\n(-x,-y,+z)\t%i\n(-x,-y,-z)\t%i\n(-x,+y,-z)\t%i\n",
                 leaving_ions[0], leaving_ions[1], leaving_ions[2], leaving_ions[3],
                 leaving_ions[4], leaving_ions[5], leaving_ions[6], leaving_ions[7]);
        write_string_to_file (str_temp, str_temp2);
        str_temp[base_name_length] = '\0';
    }

    /* sum of ints and vacs and so on for each element in each material */
    if (do_not_store_damage == 0) {
        for (i=0; i<number_of_materials; i++) {
            if (print_level >= 2) printf ("Storing for material %i:\n", i);

            if (do_not_store_damage == 0) {
                sprintf (str_temp+base_name_length, "int.mat%i.msh", i);
                if (print_level >= 2) printf (" Recoil interstitials to %s\n", str_temp);
                write_int_array_to_msh_file (str_temp, i, target_total_interstitials,
                                        list_of_materials[i].target_implanted_recoils_int, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, "repl.mat%i.msh", i);
                if (print_level >= 2) printf (" Recoil replacements to %s\n", str_temp);
                write_int_array_to_msh_file (str_temp, i, target_total_replacements,
                                        list_of_materials[i].target_implanted_recoils_repl, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, "vac.mat%i.msh", i);
                if (print_level >= 2) printf (" Vacancies to     %s\n", str_temp);
                write_int_array_to_msh_file (str_temp, i, target_total_vacancies,
                                        list_of_materials[i].target_elemental_vacancies, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, "disp.mat%i.msh", i);
                if (print_level >= 2) printf (" Displacements to   %s\n", str_temp);
                write_int_array_to_msh_file (str_temp, i, target_total_displacements,
                                        list_of_materials[i].target_elemental_disp, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';
            }

            for (j=0; j<list_of_materials[i].element_count; j++) {  /* Go through elements, store arrays for each */
                if (detailed_sputtering == 1) {
                    /* particles leaving the sample (sputtered or implanted deeper): */
                    sprintf (str_temp+base_name_length, "leaving_directions.z%i.m%.3f.mat%i.elem%i",
                             list_of_materials[i].elements_Z[j],
                             list_of_materials[i].elements_M[j], i, j);
                    if (print_level >= 2) printf (" Leaving atoms to   %s\n", str_temp);
                    sprintf (str_temp2,  "(+x,+y,+z)\t%i\n(+x,-y,+z)\t%i\n(+x,-y,-z)\t%i\n(+x,+y,-z)\t%i\n(-x,+y,+z)\t%i\n(-x,-y,+z)\t%i\n(-x,-y,-z)\t%i\n(-x,+y,-z)\t%i\n",
                             list_of_materials[i].sputter_counter[8*j+0],
                             list_of_materials[i].sputter_counter[8*j+1],
                             list_of_materials[i].sputter_counter[8*j+2],
                             list_of_materials[i].sputter_counter[8*j+3],
                             list_of_materials[i].sputter_counter[8*j+4],
                             list_of_materials[i].sputter_counter[8*j+5],
                             list_of_materials[i].sputter_counter[8*j+6],
                             list_of_materials[i].sputter_counter[8*j+7]);
                    write_string_to_file (str_temp, str_temp2);
                    str_temp[base_name_length] = '\0';
                }
            }
            if (detailed_sputtering == 1) {
                sprintf (str_temp+base_name_length, "leaving.mat%i.msh", i);
                if (print_level >= 2) printf (" Leaving from cell  %s\n", str_temp);
                write_int_array_to_msh_file (str_temp, i, target_total_sputtered,
                                    list_of_materials[i].target_sputtered_atoms, cell_count,
                                    list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';
            }
        }
    }

    if (store_transmitted_ions == 1) {
        strcat (str_temp, "transmitted.ions");
        if (print_level >= 2) printf ("Storing transmitted ions to: %s\n", str_temp);
        store_transmission_array (str_temp, transmit_list, transmission_pointer);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    }
    if (store_exiting_recoils == 1) {
        if (print_level >= 2) printf ("Storing transmitted recoils...\n");
        for (i=0; i<number_of_materials; i++) {  /* go through mats */
            for (j=0; j<list_of_materials[i].element_count; j++) { /* go through elements, store arrays for each */
                sprintf (str_temp+base_name_length, "leaving_recoils.z%i.m%.3f.mat%i.elem%i",
                         list_of_materials[i].elements_Z[j],
                         list_of_materials[i].elements_M[j], i, j);
                store_transmission_array (str_temp, list_of_materials[i].elemental_leaving_recoils[j],
                                          list_of_materials[i].leaving_recoils_pointer[j]);
                str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
            }
        }
    }

    /* depth distribution functions */
    //if (store_depth_dist == 1) {
        do_depth_dist_statistics ();
        strcat (str_temp, "depth_dist_functions.dat");
        store_depth_dist_array (str_temp);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */

        do_radial_dist_statistics ();
        strcat (str_temp, "radial_dist_functions.dat");
        store_radial_dist_array (str_temp);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    //}

    free (str_temp);

    return 0;
}

/*=============================================================================
  Function Name : write_int_array_to_msh_file
  Description   : Writes the designated array of Count elements into a msh file.

  Inputs  :
            char* file_name
            int mat_i
            int* source_array_total
            int *source_array_part
            int count
            int n_element
  Outputs : no.

  Notes :
      The caller needs to make sure, that the array is large enough.

      Function call:
        target - int get_target_XYZ (int index, int* x, int* y, int* z);
=============================================================================*/
int write_int_array_to_msh_file (char *file_name, int mat_i, int *source_array_total,
                                 int *source_array_part[], int count, int n_element) {
    FILE *fp;
    int i, j, k, jj;
    int num_nodes, num_elements;
    int  x, y, z;
    int  *element_non_zero, * element_non_zero_index, *node_non_zero, *node_non_zero_index;
    float rx, ry, rz;
    double node_value;

    /* only count and output the non-zero elements */
    element_non_zero = (int*) calloc (count, sizeof (int));
    if (element_non_zero == NULL) return -6001;  /* cannot allocate memory */
    element_non_zero_index = (int*) calloc (count, sizeof (int));
    if (element_non_zero_index == NULL) return -6002;  /* cannot allocate memory */
    for (i=0; i<count; i++) {
        element_non_zero[i] = -1;
    }
    num_elements = 0;
    for (i=0; i<count; i++) {
        if (source_array_total[i] > 0.0) {
            element_non_zero_index[num_elements] = i;
            num_elements ++;
            element_non_zero[i] = num_elements;
        }
        else {
            for (j=0; j<n_element; j++) {
                if (source_array_part [j][i] > 0.0) {
                    element_non_zero_index[num_elements] = i;
                    num_elements ++;
                    element_non_zero[i] = num_elements;
                    break;
                }
            }
        }
    }

    /* count notes */
    node_non_zero = (int*) calloc (node_count, sizeof (int));
    if (node_non_zero == NULL) return -6003;  /* cannot allocate memory */
    node_non_zero_index = (int*) calloc (node_count, sizeof (int));
    if (node_non_zero_index == NULL) return -6004;  /* cannot allocate memory */
    num_nodes = 0;
    for (i=0; i<node_count; i++) {
        get_node_XYZ (i, &x, &y, &z);

        /* (x, y, z) */
        k = get_target_index (x, y, z);
        if (element_non_zero[k] >=0) goto out5000;
        /* (x-1, y, z) */
        if (x-1 >= 0) {
            k = get_target_index (x-1, y, z);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x, y-1, z) */
        if (y-1 >= 0) {
            k = get_target_index (x, y-1, z);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x, y, z-1) */
        if (z-1 >= 0) {
            k = get_target_index (x, y, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x-1, y-1, z) */
        if (x-1>= 0 && y-1>= 0) {
            k = get_target_index (x-1, y-1, z);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x-1, y, z-1) */
        if (x-1>= 0 && z-1>= 0) {
            k = get_target_index (x-1, y, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x, y-1, z-1) */
        if (y-1>= 0 && z-1>= 0) {
            k = get_target_index (x, y-1, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x-1, y-1, z-1) */
        if (x-1>=0 && y-1>=0 && z-1>=0) {
            k = get_target_index (x-1, y-1, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }

        continue;

out5000:
        node_non_zero_index[num_nodes] = i;
        num_nodes ++;
        node_non_zero[i] = num_nodes;
    }

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -6005;

    /* MeshFormat */
    fprintf (fp, "$MeshFormat\n");
    //if (msh_file_type = 0) {  /* ASCII */
        fprintf (fp, "2.2 0 8\n");  /* version-number, file-type, data-size */
    //}
    //else {  /* TODO: binary */

    //}
    fprintf (fp, "$EndMeshFormat\n");

    /* Nodes */
    fprintf (fp, "$Nodes\n");
    fprintf (fp, "%i\n", num_nodes);  /* number-of-nodes */
    for (i=0; i<num_nodes; i++) {
        get_node_XYZ (node_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        get_node_relative_XYZ (x, y, z, &rx, &ry, &rz);

        /* node-number, x-coord, y-coord, z-coord */
        fprintf (fp, "%i\t%g\t%g\t%g\n", i+1, rx, ry, rz);
    }
    fprintf (fp, "$EndNodes\n");

    /* Elements */
    fprintf (fp, "$Elements\n");
    fprintf (fp, "%i\n", num_elements);
    for (i=0; i<num_elements; i++) {
        get_target_XYZ (element_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        fprintf (fp, "%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", i+1, 5, 1, 1,
                  node_non_zero[get_node_index (x, y, z)],       node_non_zero[get_node_index (x, y+1, z)]    ,
                  node_non_zero[get_node_index (x, y+1, z+1)],   node_non_zero[get_node_index (x, y, z+1)]    ,
                  node_non_zero[get_node_index (x+1, y, z)],     node_non_zero[get_node_index (x+1, y+1, z)]  ,
                  node_non_zero[get_node_index (x+1, y+1, z+1)], node_non_zero[get_node_index (x+1, y, z+1)] );
    }
    fprintf (fp, "$EndElements\n");

    /* ElementData */
    fprintf (fp, "$NodeData\n");
    fprintf (fp, "%i\n", 1);  /* one string tag */
    fprintf (fp, "Mat%i_Total_node\n", mat_i);  /* the name of the view */
    fprintf (fp, "%i\n", 1);  /* one real tag */
    fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    fprintf (fp, "%i\n", 3);  /* three integer tags */
    fprintf (fp, "%i\n", 0);  /* the time step */
    fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    fprintf (fp, "%i\n", num_nodes);  /* num_elements associated elemental values */
    for (i=0; i<num_nodes; i++) {
        get_node_XYZ (node_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        /* (x, y, z) */
        k = 0;
        node_value = 0.0;
        j = get_target_index (x, y, z);
        if (element_non_zero[j] >= 0) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x-1, y, z) */
        j = get_target_index (x-1, y, z);
        if (x-1 >= 0 && element_non_zero[j] >= 0) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x, y-1, z) */
        j = get_target_index (x, y-1, z);
        if (y-1 >= 0 && element_non_zero[j] >=0) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x, y, z-1) */
        j = get_target_index (x, y, z-1);
        if (z-1 >= 0 && element_non_zero[j] >= 0) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x-1, y-1, z) */
        j = get_target_index (x-1, y-1, z);
        if (x-1>= 0 && y-1>= 0 && element_non_zero[j >= 0]) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x-1, y, z-1) */
        j = get_target_index (x-1, y, z-1);
        if (x-1>= 0 && z-1>= 0 && element_non_zero[j] >= 0) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x, y-1, z-1) */
        j = get_target_index (x, y-1, z-1);
        if (y-1>= 0 && z-1>= 0 && element_non_zero[j] >=0) {
            node_value += source_array_total[j];
            k ++;
        }
        /* (x-1, y-1, z-1) */
        j = get_target_index (x-1, y-1, z-1);
        if (x-1>=0 && y-1>=0 && z-1>=0 && element_non_zero[j] >=0) {
            node_value += source_array_total[j];
            k ++;
        }
        node_value /= k;
        if (normalize_output == 1)
            fprintf (fp, "%i\t%g\n", i+1, (double) node_value * unit_conversion_factor);
        else {
            fprintf (fp, "%i\t%g\n", i+1, (double) node_value);
        }
    }

    for (jj=0; jj<n_element; jj++) {
        fprintf (fp, "$NodeData\n");
        fprintf (fp, "%i\n", 1);  /* one string tag */
        fprintf (fp, "\"Mat%i_Element-%i_node\"\n", mat_i, jj+1);  /* the name of the view */
        fprintf (fp, "%i\n", 1);  /* one real tag */
        fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
        fprintf (fp, "%i\n", 3);  /* three integer tags */
        fprintf (fp, "%i\n", 0);  /*the time step */
        fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
        fprintf (fp, "%i\n", num_nodes);  /* num_elements associated elemental values */
        for (i=0; i<num_nodes; i++) {
            get_node_XYZ (node_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
            /* (x, y, z) */
            k = 0;
            node_value = 0.0;
            j = get_target_index (x, y, z);
            if (element_non_zero[j] >= 0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x-1, y, z) */
            j = get_target_index (x-1, y, z);
            if (x-1 >= 0 && element_non_zero[j] >= 0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x, y-1, z) */
            j = get_target_index (x, y-1, z);
            if (y-1 >= 0 && element_non_zero[j] >=0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x, y, z-1) */
            j = get_target_index (x, y, z-1);
            if (z-1 >= 0 && element_non_zero[j] >= 0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x-1, y-1, z) */
            j = get_target_index (x-1, y-1, z);
            if (x-1>= 0 && y-1>= 0 && element_non_zero[j >= 0]) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x-1, y, z-1) */
            j = get_target_index (x-1, y, z-1);
            if (x-1>= 0 && z-1>= 0 && element_non_zero[j] >= 0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x, y-1, z-1) */
            j = get_target_index (x, y-1, z-1);
            if (y-1>= 0 && z-1>= 0 && element_non_zero[j] >=0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            /* (x-1, y-1, z-1) */
            j = get_target_index (x-1, y-1, z-1);
            if (x-1>=0 && y-1>=0 && z-1>=0 && element_non_zero[j] >=0) {
                node_value += source_array_part[jj][j];
                k ++;
            }
            node_value /= k;
            if (normalize_output == 1)
                fprintf (fp, "%i\t%g\n", i+1, (double) node_value * unit_conversion_factor);
            else {
                fprintf (fp, "%i\t%g\n", i+1, (double) node_value);
            }
        }
    }

    /* ElementData */
    //fprintf (fp, "$ElementData\n");
    //fprintf (fp, "%i\n", 1);  /* one string tag */
    //fprintf (fp, "Mat%i_Total_element\n", mat_i);  /* the name of the view */
    //fprintf (fp, "%i\n", 1);  /* one real tag */
    //fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    //fprintf (fp, "%i\n", 3);  /* three integer tags */
    //fprintf (fp, "%i\n", 0);  /* the time step */
    //fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    //fprintf (fp, "%i\n", num_elements);  /* num_elements associated elemental values */
    //if (normalize_output == 1) {
    //    for (i=0; i<num_elements; i++)
    //        fprintf (fp, "%i\t%g\n", i+1, (double) source_array_total[i] * unit_conversion_factor);
    //}
    //else {
    //    for (i=0; i<num_elements; i++)
    //        fprintf (fp, "%i\t%g\n", i+1, (double) source_array_total[i]);
    //}

    //for (j=0; j<n_element; j++) {
    //    fprintf (fp, "$ElementData\n");
    //    fprintf (fp, "%i\n", 1);  /* one string tag */
    //    fprintf (fp, "\"Mat%i_Element-%i_element\"\n", mat_i, j+1);  /* the name of the view */
    //    fprintf (fp, "%i\n", 1);  /* one real tag */
    //    fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    //    fprintf (fp, "%i\n", 3);  /* three integer tags */
    //    fprintf (fp, "%i\n", 0);  /*the time step */
    //    fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    //    fprintf (fp, "%i\n", num_elements);  /* num_elements associated elemental values */
    //    if (normalize_output == 1) {
    //        for (i=0; i<num_elements; i++)
    //            fprintf (fp, "%i\t%g\n", i+1, (double) source_array_part[j][i] * unit_conversion_factor);
    //    }
    //    else {
    //        for (i=0; i<num_elements; i++)
    //            fprintf (fp, "%i\t%g\n", i+1, (double) source_array_part[j][i]);
    //    }
    //}

    fclose (fp);

    free (element_non_zero);
    free (element_non_zero_index);
    free (node_non_zero);
    free (node_non_zero_index);

    return 0;
}

/*=============================================================================
  Function Name : write_double_array_to_msh_file
  Description   : Writes the designated array of Count elements into a msh file.

  Inputs  :
            char* file_name
            int mat_i
            int* source_array_total
            int *source_array_part
            int count
            int n_element
  Outputs : no.

  Notes :
      The caller needs to make sure, that the array is large enough.

      Function call:
        target - int get_target_XYZ (int index, int* x, int* y, int* z);
=============================================================================*/
int write_double_array_to_msh_file (char *file_name, double *source_array1,
                                    double *source_array2, int count) {
    FILE *fp;
    int i, j, k;
    int num_nodes, num_elements;
    int  x, y, z;
    int  *element_non_zero, *element_non_zero_index, *node_non_zero, *node_non_zero_index;
    float rx, ry, rz;
    double node_value;

    /* only count and output the non-zero elements */
    element_non_zero = (int*) calloc (count, sizeof (int));
    if (element_non_zero == NULL) return -6006;  /* cannot allocate memory */
    element_non_zero_index = (int*) calloc (count, sizeof (int));
    if (element_non_zero_index == NULL) return -6007;  /* cannot allocate memory */
    for (i=0; i<count; i++) {
        element_non_zero[i] = -1;
    }
    num_elements = 0;
    for (i=0; i<count; i++) {
        if (source_array1[i] > 0.0 || source_array2[i] > 0.0) {
            element_non_zero_index[num_elements] = i;
            num_elements ++;
            element_non_zero[i] = num_elements;
        }
    }

    /* count notes */
    node_non_zero = (int*) calloc (node_count, sizeof (int));
    if (node_non_zero == NULL) return -6008;  /* cannot allocate memory */
    node_non_zero_index = (int*) calloc (node_count, sizeof (int));
    if (node_non_zero_index == NULL) return -6009;  /* cannot allocate memory */
    num_nodes = 0;
    for (i=0; i<node_count; i++) {
        get_node_XYZ (i, &x, &y, &z);

        /* (x, y, z) */
        k = get_target_index (x, y, z);
        if (element_non_zero[k] >=0) goto out5000;
        /* (x-1, y, z) */
        if (x-1 >= 0) {
            k = get_target_index (x-1, y, z);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x, y-1, z) */
        if (y-1 >= 0) {
            k = get_target_index (x, y-1, z);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x, y, z-1) */
        if (z-1 >= 0) {
            k = get_target_index (x, y, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x-1, y-1, z) */
        if (x-1>= 0 && y-1>= 0) {
            k = get_target_index (x-1, y-1, z);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x-1, y, z-1) */
        if (x-1>= 0 && z-1>= 0) {
            k = get_target_index (x-1, y, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x, y-1, z-1) */
        if (y-1>= 0 && z-1>= 0) {
            k = get_target_index (x, y-1, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }
        /* (x-1, y-1, z-1) */
        if (x-1>=0 && y-1>=0 && z-1>=0) {
            k = get_target_index (x-1, y-1, z-1);
            if (element_non_zero[k] >=0) goto out5000;
        }

        continue;

out5000:
        node_non_zero_index[num_nodes] = i;
        num_nodes ++;
        node_non_zero[i] = num_nodes;
    }

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -6010;

    /* MeshFormat */
    fprintf (fp, "$MeshFormat\n");
    //if (file_type = 0) {  /* ASCII */
        fprintf (fp, "2.2 0 8\n");  /* version-number, file-type, data-size */
    //}
    //else {  /* TODO: binary */

    //}
    fprintf (fp, "$EndMeshFormat\n");

    /* Nodes */
    fprintf (fp, "$Nodes\n");
    fprintf (fp, "%i\n", num_nodes);  /* number-of-nodes */
    for (i=0; i<num_nodes; i++) {
        get_node_XYZ (node_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        get_node_relative_XYZ (x, y, z, &rx, &ry, &rz);

        /* node-number, x-coord, y-coord, z-coord */
        fprintf (fp, "%i\t%g\t%g\t%g\n", i+1, rx, ry, rz);
    }
    fprintf (fp, "$EndNodes\n");

    /* Elements */
    fprintf (fp, "$Elements\n");
    fprintf (fp, "%i\n", num_elements);
    for (i=0; i<num_elements; i++) {
        get_target_XYZ (element_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        fprintf (fp, "%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", i+1, 5, 1, 1,
                 node_non_zero[get_node_index (x, y, z)],       node_non_zero[get_node_index (x, y+1, z)]    ,
                 node_non_zero[get_node_index (x, y+1, z+1)],   node_non_zero[get_node_index (x, y, z+1)]    ,
                 node_non_zero[get_node_index (x+1, y, z)],     node_non_zero[get_node_index (x+1, y+1, z)]  ,
                 node_non_zero[get_node_index (x+1, y+1, z+1)], node_non_zero[get_node_index (x+1, y, z+1)] );
                 //if (node_non_zero[get_node_index (x, y, z)] == 0) printf ("%i\t%i\t%i\t%i\n", get_node_index (x, y, z), x, y, z);
    }
    fprintf (fp, "$EndElements\n");

    /* NodeData */
    fprintf (fp, "$NodeData\n");
    fprintf (fp, "%i\n", 1);  /* one string tag */
    fprintf (fp, "\"target_energy_electrons_node\"\n");  /* the name of the view */
    fprintf (fp, "%i\n", 1);  /* one real tag */
    fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    fprintf (fp, "%i\n", 3);  /* three integer tags */
    fprintf (fp, "%i\n", 0);  /* the time step */
    fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    fprintf (fp, "%i\n", num_nodes);  /* num_elements associated elemental values */
    for (i=0; i<num_nodes; i++) {
        get_node_XYZ (node_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        /* (x, y, z) */
        k = 0;
        node_value = 0.0;
        j = get_target_index (x, y, z);
        if (element_non_zero[j] >= 0) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x-1, y, z) */
        j = get_target_index (x-1, y, z);
        if (x-1 >= 0 && element_non_zero[j] >= 0) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x, y-1, z) */
        j = get_target_index (x, y-1, z);
        if (y-1 >= 0 && element_non_zero[j] >=0) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x, y, z-1) */
        j = get_target_index (x, y, z-1);
        if (z-1 >= 0 && element_non_zero[j] >= 0) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x-1, y-1, z) */
        j = get_target_index (x-1, y-1, z);
        if (x-1>= 0 && y-1>= 0 && element_non_zero[j >= 0]) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x-1, y, z-1) */
        j = get_target_index (x-1, y, z-1);
        if (x-1>= 0 && z-1>= 0 && element_non_zero[j] >= 0) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x, y-1, z-1) */
        j = get_target_index (x, y-1, z-1);
        if (y-1>= 0 && z-1>= 0 && element_non_zero[j] >=0) {
            node_value += source_array1[j];
            k ++;
        }
        /* (x-1, y-1, z-1) */
        j = get_target_index (x-1, y-1, z-1);
        if (x-1>=0 && y-1>=0 && z-1>=0 && element_non_zero[j] >=0) {
            node_value += source_array1[j];
            k ++;
        }
        node_value /= k;
        if (normalize_output == 1)
            fprintf (fp, "%i\t%g\n", i+1, (double) node_value * unit_conversion_factor);
        else {
            fprintf (fp, "%i\t%g\n", i+1, (double) node_value);
        }
    }

    /* NodeData */
    fprintf (fp, "$NodeData\n");
    fprintf (fp, "%i\n", 1);  /* one string tag */
    fprintf (fp, "\"target_energy_phonons_node\"\n");  /* the name of the view */
    fprintf (fp, "%i\n", 1);  /* one real tag */
    fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    fprintf (fp, "%i\n", 3);  /* three integer tags */
    fprintf (fp, "%i\n", 0);  /* the time step */
    fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    fprintf (fp, "%i\n", num_nodes);  /* num_elements associated elemental values */
    for (i=0; i<num_nodes; i++) {
        get_node_XYZ (node_non_zero_index[i], &x, &y, &z);  /* convert linear index to 3 coords */
        /* (x, y, z) */
        k = 0;
        node_value = 0.0;
        j = get_target_index (x, y, z);
        if (element_non_zero[j] >= 0) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x-1, y, z) */
        j = get_target_index (x-1, y, z);
        if (x-1 >= 0 && element_non_zero[j] >= 0) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x, y-1, z) */
        j = get_target_index (x, y-1, z);
        if (y-1 >= 0 && element_non_zero[j] >=0) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x, y, z-1) */
        j = get_target_index (x, y, z-1);
        if (z-1 >= 0 && element_non_zero[j] >= 0) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x-1, y-1, z) */
        j = get_target_index (x-1, y-1, z);
        if (x-1>= 0 && y-1>= 0 && element_non_zero[j >= 0]) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x-1, y, z-1) */
        j = get_target_index (x-1, y, z-1);
        if (x-1>= 0 && z-1>= 0 && element_non_zero[j] >= 0) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x, y-1, z-1) */
        j = get_target_index (x, y-1, z-1);
        if (y-1>= 0 && z-1>= 0 && element_non_zero[j] >=0) {
            node_value += source_array2[j];
            k ++;
        }
        /* (x-1, y-1, z-1) */
        j = get_target_index (x-1, y-1, z-1);
        if (x-1>=0 && y-1>=0 && z-1>=0 && element_non_zero[j] >=0) {
            node_value += source_array2[j];
            k ++;
        }
        node_value /= k;
        if (normalize_output == 1) {
            fprintf (fp, "%i\t%g\n", i+1, (double) node_value * unit_conversion_factor);
        }
        else {
            fprintf (fp, "%i\t%g\n", i+1, (double) node_value);
        }
    }

    /* ElementData */
    //fprintf (fp, "$ElementData\n");
    //fprintf (fp, "%i\n", 1);  /* one string tag */
    //fprintf (fp, "\"target_energy_electrons_element\"\n");  /* the name of the view */
    //fprintf (fp, "%i\n", 1);  /* one real tag */
    //fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    //fprintf (fp, "%i\n", 3);  /* three integer tags */
    //fprintf (fp, "%i\n", 0);  /* the time step */
    //fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    //fprintf (fp, "%i\n", num_elements);  /* num_elements associated elemental values */
    //if (normalize_output == 1) {
    //    for (i=0; i<num_elements; i++)
    //        fprintf (fp, "%i\t%g\n", i+1, (double) source_array1[i] * unit_conversion_factor);
    //}
    //else {
    //    for (i=0; i<num_elements; i++)
    //        fprintf (fp, "%i\t%g\n", i+1, (double) source_array1[i]);
    //}

    //fprintf (fp, "$ElementData\n");
    //fprintf (fp, "%i\n", 1);  /* one string tag */
    //fprintf (fp, "\"target_energy_phonons_element\"\n");  /* the name of the view */
    //fprintf (fp, "%i\n", 1);  /* one real tag */
    //fprintf (fp, "%g\n", (double) max_no_ions);  /* the time value */
    //fprintf (fp, "%i\n", 3);  /* three integer tags */
    //fprintf (fp, "%i\n", 0);  /*the time step */
    //fprintf (fp, "%i\n", 1);  /* 1-component (scalar) field */
    //fprintf (fp, "%i\n", num_elements);  /* num_elements associated elemental values */
    //if (normalize_output == 1) {
    //    for (i=0; i<num_elements; i++)
    //        fprintf (fp, "%i\t%g\n", i+1, (double) source_array2[i] * unit_conversion_factor);
    //}
    //else {
    //    for (i=0; i<num_elements; i++)
    //        fprintf (fp, "%i\t%g\n", i+1, (double) source_array2[i]);
    //}

    fclose (fp);

    free (element_non_zero);
    free (element_non_zero_index);
    free (node_non_zero);
    free (node_non_zero_index);

    return 0;
}

/*=============================================================================
  Function Name : get_node_index
  Description   : Calculates index for accessing one-dimensional node arrays
                  knowing the integer cell coordinates.

  Inputs  : int x, int y, int z
  Outputs : no.

  Notes :
         layer_count_xy = cell_count_x * cell_count_y,
         still use layer_count_yz is in order to make output in style of
         (x, y, z, value).
=============================================================================*/
int get_node_index (int x, int y, int z) {
    return x * node_yz + y * node_z + z;
}

/*=============================================================================
  Function Name : get_node_XYZ
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
void get_node_XYZ (int index, int *x, int *y, int *z) {
    *z = ( index % node_yz) % node_z;
    *y = ((index % node_yz) - (*z)) / node_z;
    *x = ((index - (*z) - node_z * (*y)) / node_yz);

    return;
}

/*=============================================================================
  Function Name : get_node_relative_XYZ
  Description   : Calculates the relative x,y,z coordinates (unit: traget_size)
                  for a given position.

  Inputs  : int x, int y, int z

  Outputs : float *rx, float *ry, float *rz

  Notes : no.
=============================================================================*/
void get_node_relative_XYZ (int x, int y, int z, float *rx, float *ry, float *rz) {

    *rx = (float) x / (float) node_x;
    *ry = (float) y / (float) node_y;
    *rz = (float) z / (float) node_z;

    return;
}
