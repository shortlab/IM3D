/******************************************************************************
  Module Name : cfgwriter.c
  Module Date : 09/03/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Output files with .cfg file format.

  Others :
      Error numbers in this module 5000-5999.

      Revision History:
      Date    Rel Ver.    Notes
      04/01/2014  1.0.0  first coded and integrated in module fileio
      09/03/2014  1.0.4  rewritten in module cfgwriter
******************************************************************************/
#include "cfgwriter.h"

/*=============================================================================
  Function Name : store_results_cfg
  Description   : Store the results of the simulation (arrays with distribution
                  of implanted ions, defects etc.) in .cfg format.
                  This is the non-dynamic version of the function for
                  material-based output.

  Inputs  :
            char *base_name
  Outputs : no.

  Notes : ygli
=============================================================================*/
int store_results_cfg (char *base_name) {
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
    strcat (str_temp, "ions.total.cfg");
    if (print_level >= 2) printf ("Storing implanted ions to:      %s\n", str_temp);
    write_int_array_to_cfg_file (str_temp, 0, target_implanted_ions, 0, cell_count, 0);
    str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */

    if (do_not_store_damage == 0) {
        /* Part of ions that replaced identical target atoms */
        strcat (str_temp, "ions.replacements.cfg");
        if (print_level >= 2) printf ("Storing implanted replacing ions to:   %s\n", str_temp);
        write_int_array_to_cfg_file (str_temp, 0, target_replacing_ions, 0, cell_count, 0);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    }

    sum_up_material_arrays ();  /* sum up the individual element-dependent arrays */

    /* deposited energy */
    if (store_energy_deposit == 1) {
        strcat (str_temp, "energy.deposit.cfg");
        if (print_level >= 2) printf ("Storing energy deposited to electrons and phonons:  %s\n", str_temp);
        write_double_array_to_cfg_file (str_temp, target_energy_electrons, target_energy_phonons, cell_count);
        str_temp[base_name_length] = '\0';
    }

    /* store single ion sputter yields */
    if (single_ion_sputter_yields == 1) {
        strcat (str_temp, "single_ion_sputter_yields");
        if (print_level >= 2) printf ("Storing single ion sputter yields to:   %s\n", str_temp);
        fp = fopen (str_temp, "w");
        if (fp == NULL) {
            printf("Error: cannot store histogram.\n");
            return -5000;
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
                sprintf (str_temp+base_name_length, "int.mat%i.cfg", i);
                if (print_level >= 2) printf (" Recoil interstitials to %s\n", str_temp);
                write_int_array_to_cfg_file (str_temp, i, target_total_interstitials,
                                        list_of_materials[i].target_implanted_recoils_int, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, "repl.mat%i.cfg", i);
                if (print_level >= 2) printf (" Recoil replacements to %s\n", str_temp);
                write_int_array_to_cfg_file (str_temp, i, target_total_replacements,
                                        list_of_materials[i].target_implanted_recoils_repl, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, "vac.mat%i.cfg", i);
                if (print_level >= 2) printf (" Vacancies to     %s\n", str_temp);
                write_int_array_to_cfg_file (str_temp, i, target_total_vacancies,
                                        list_of_materials[i].target_elemental_vacancies, cell_count,
                                        list_of_materials[i].element_count);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, "disp.mat%i.cfg", i);
                if (print_level >= 2) printf (" Displacements to   %s\n", str_temp);
                write_int_array_to_cfg_file (str_temp, i, target_total_displacements,
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
                sprintf (str_temp+base_name_length, "leaving.mat%i.cfg", i);
                if (print_level >= 2) printf (" Leaving from cell  %s\n", str_temp);
                write_int_array_to_cfg_file (str_temp, i,  target_total_sputtered,
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
  Function Name : write_int_array_to_cfg_file
  Description   : Writes the designated array of Count elements into a cfg file.
                  If the file is just one column of values then FileType should
                  be set to 1. If the file contains 4 columns (x, y, z, value)
                  then set it to 0.

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
int write_int_array_to_cfg_file (char *file_name, int mat_i, int *source_array_total,
                                 int *source_array_part[], int count, int n_element) {
    FILE *fp;
    int  i, j, count_non_zero;
    int  x, y, z;
    int  *flag;
    float rx, ry, rz;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -5001;

    /* only count and output the non-zero elements */
    flag = (int*) calloc (count, sizeof (int));
    if (flag == NULL) return -5002;  /* cannot allocate memory */
    count_non_zero = 0;
    for (i=0; i<count; i++) {
        if (source_array_total[i] > 0) {
            count_non_zero ++;
            flag[i] = 1;
        }
        else {
            for (j=0; j<n_element; j++) {
                if (source_array_part [j][i] > 0) {
                    count_non_zero ++;
                    flag[i] = 1;
                    break;
                }
            }
        }
    }

    /* .cfg format */
    //fprintf (fp, "Number of particles = %i\n", count);
    fprintf (fp, "Number of particles = %i\n", count_non_zero);
    //fprintf (fp, "A = 1 nm (basic length-scale)\n");
    fprintf (fp, "A = 10 Angstrom (basic length-scale)\n");
    fprintf (fp, "H0(1,1) = %g A\n", target_size_x);
    fprintf (fp, "H0(1,2) = 0.0 A\n");
    fprintf (fp, "H0(1,3) = 0.0 A\n");
    fprintf (fp, "H0(2,1) = 0.0 A\n");
    fprintf (fp, "H0(2,2) = %g A\n", target_size_y);
    fprintf (fp, "H0(2,3) = 0.0 A\n");
    fprintf (fp, "H0(3,1) = 0.0 A\n");
    fprintf (fp, "H0(3,2) = 0.0 A\n");
    fprintf (fp, "H0(3,3) = %g A\n", target_size_z);
    fprintf (fp, ".NO_VELOCITY.\n");
    fprintf (fp, "entry_count = %i\n", n_element + 4);
    fprintf (fp, "auxiliary[0] = Mat%i_total\n", mat_i);
    if (n_element >= 1) {
        fprintf (fp, "auxiliary[1] = Z%i\n", list_of_materials[mat_i].elements_Z[0]);
    }
    if (n_element >= 2) {
        fprintf (fp, "auxiliary[2] = Z%i\n", list_of_materials[mat_i].elements_Z[1]);
    }
    if (n_element >= 3) {
        fprintf (fp, "auxiliary[3] = Z%i\n", list_of_materials[mat_i].elements_Z[2]);
    }
    if (n_element >= 4) {
        fprintf (fp, "auxiliary[4] = Z%i\n", list_of_materials[mat_i].elements_Z[3]);
    }
    if (n_element >= 5) {
        fprintf (fp, "auxiliary[5] = Z%i\n", list_of_materials[mat_i].elements_Z[4]);
    }
    if (n_element >= 6) {
        fprintf (fp, "ERROR: Print too many elements in cfg file, no more than 5 elements!\n");
    }
    fprintf (fp, "%i\n", n_element);
    fprintf (fp, "%i\n", mat_i);

    /* .cfg form, x,y,z and value. */
    if (normalize_output == 1) {
        for (i=0; i<count; i++) {
            /* only count and output the non-zero elements */
            if (flag[i] == 1) {
                get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
                get_relative_XYZ (x, y, z, &rx, &ry, &rz);
                switch (n_element) {
                case 0 :
                    fprintf (fp, "%f\t%f\t%f\t%f\n", rx, ry, rz, source_array_total[i] * unit_conversion_factor);
                    break;
                case 1 :
                    fprintf (fp, "%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array_total[i] * unit_conversion_factor,
                        source_array_part [i][0] * unit_conversion_factor);
                    break;
                case 2 :
                    fprintf (fp, "%f\t%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array_total[i] * unit_conversion_factor,
                        source_array_part [0][i] * unit_conversion_factor, source_array_part [1][i] * unit_conversion_factor);
                    break;
                case 3 :
                    fprintf (fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array_total[i] * unit_conversion_factor,
                        source_array_part [0][i] * unit_conversion_factor, source_array_part [1][i] * unit_conversion_factor,
                        source_array_part [2][i] * unit_conversion_factor);
                    break;
                case 4 :
                    fprintf (fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array_total[i] * unit_conversion_factor,
                        source_array_part [0][i] * unit_conversion_factor, source_array_part [1][i] * unit_conversion_factor,
                        source_array_part [2][i] * unit_conversion_factor, source_array_part [3][i] * unit_conversion_factor);
                    break;
                case 5 :
                    fprintf (fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array_total[i] * unit_conversion_factor,
                        source_array_part [0][i] * unit_conversion_factor, source_array_part [1][i] * unit_conversion_factor,
                        source_array_part [2][i] * unit_conversion_factor, source_array_part [3][i] * unit_conversion_factor,
                        source_array_part [4][i] * unit_conversion_factor);
                    break;
                default:
                    printf ("ERROR1: Print too many elements in cfg file, no more than 5 elements!\n");
                    break;
                }
            }
        }
    }
    else {
        for (i=0; i<count; i++) {
            /* only count and output the non-zero elements */
            if (flag[i] == 1) {
                get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
                get_relative_XYZ (x, y, z, &rx, &ry, &rz);
                switch (n_element) {
                case 0 :
                    fprintf (fp, "%f\t%f\t%f\t%i\n", rx, ry, rz, source_array_total[i]);
                    break;
                case 1 :
                    fprintf (fp, "%f\t%f\t%f\t%i\t%i\n", rx, ry, rz, source_array_total[i], source_array_part [0][i]);
                    break;
                case 2 :
                    fprintf (fp, "%f\t%f\t%f\t%8i\t%8i\t%8i\n", rx, ry, rz, source_array_total[i],
                             source_array_part [0][i], source_array_part [1][i]);
                    break;
                case 3 :
                    fprintf (fp, "%f\t%f\t%f\t%i\t%i\t%i\t%i\n", rx, ry, rz, source_array_total[i],
                             source_array_part [0][i], source_array_part [1][i],source_array_part [2][i]);
                    break;
                case 4 :
                    fprintf (fp, "%f\t%f\t%f\t%i\t%i\t%i\t%i\t%i\n", rx, ry, rz, source_array_total[i],
                             source_array_part [0][i], source_array_part [1][i], source_array_part [2][i],
                             source_array_part [3][i]);
                    break;
                case 5 :
                    fprintf (fp, "%f\t%f\t%f\t%i\t%i\t%i\t%i\t%i\t%i\n", rx, ry, rz, source_array_total[i],
                             source_array_part [0][i], source_array_part [1][i], source_array_part [2][i],
                             source_array_part [3][i], source_array_part [4][i]);
                    break;
                default:
                    printf ("ERROR2: Print too many elements in cfg file, no more than 5 elements!\n");
                    break;
                }
            }
        }
    }

    fclose (fp);

    free (flag);

    return 0;
}

/*=============================================================================
  Function Name : write_double_array_to_cfg_file
  Description   : Writes the double array of energy deposited into a .cfg file.
                  If the file is just one column of values then FileType should
                  be set to 1. If the file contains 4 columns (x, y, z, value)
                  then set it to 0.

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
int write_double_array_to_cfg_file (char *file_name, double *source_array1,
                                    double *source_array2, int count) {
    FILE *fp;
    int  i, count_non_zero;
    int  x, y, z;
    float rx, ry, rz;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -5003;

    /* only count and output the non-zero elements */
    count_non_zero = 0;
    for (i=0; i<count; i++)
        if (source_array1[i]>0 || source_array2[i]>0) count_non_zero ++;

    /* .cfg format */
    //fprintf (fp, "Number of particles = %i\n", count);
    fprintf (fp, "Number of particles = %i\n", count_non_zero);
    //fprintf (fp, "A = 1 nm (basic length-scale)\n");
    fprintf (fp, "A = 10 Angstrom (basic length-scale)\n");
    fprintf (fp, "H0(1,1) = %g A\n", target_size_x);
    fprintf (fp, "H0(1,2) = 0.0 A\n");
    fprintf (fp, "H0(1,3) = 0.0 A\n");
    fprintf (fp, "H0(2,1) = 0.0 A\n");
    fprintf (fp, "H0(2,2) = %g A\n", target_size_y);
    fprintf (fp, "H0(2,3) = 0.0 A\n");
    fprintf (fp, "H0(3,1) = 0.0 A\n");
    fprintf (fp, "H0(3,2) = 0.0 A\n");
    fprintf (fp, "H0(3,3) = %g A\n", target_size_z);
    fprintf (fp, ".NO_VELOCITY.\n");
    fprintf (fp, "entry_count = %i\n", 5);
    fprintf (fp, "auxiliary[0] = electrons (eV)\n");
    fprintf (fp, "auxiliary[1] = phonons (eV)\n");
    fprintf (fp, "%i\n", 2);
    fprintf (fp, "energy_deposit\n");

    /* .cfg form, x,y,z and value. */
    if (normalize_output == 1) {
        for (i=0; i<count; i++) {
            /* only count and output the non-zero elements */
            if (source_array1[i]>1.0e-15 || source_array2[i]>1.0e-15) {
                get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
                get_relative_XYZ (x, y, z, &rx, &ry, &rz);

                fprintf (fp, "%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array1[i] * unit_conversion_factor,
                         source_array2[i] * unit_conversion_factor);
            }
        }
    }
    else {
        for (i=0; i<count; i++) {
            /* only count and output the non-zero elements */
            if (source_array1[i]>0 || source_array2[i]>0) {
                get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
                get_relative_XYZ (x, y, z, &rx, &ry, &rz);

                fprintf (fp, "%f\t%f\t%f\t%f\t%f\n", rx, ry, rz, source_array1[i], source_array2[i]);
            }
        }
    }

    fclose (fp);

    return 0;
}
