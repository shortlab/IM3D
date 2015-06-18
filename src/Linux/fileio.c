/******************************************************************************
  Module Name : fileio.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Contains the some helping function for file access, etc.

  Others :
      Error numbers in this module 4000-4999.
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
      04/01/14            add subroutine store_results_cfg;
******************************************************************************/
#include "fileio.h"

/*=============================================================================
  Function Name : read_init_file
  Description   : A general reader for the various input config files. It parses
                  and ini-like file. Whenever it finds a new DataBlock it calls
                  the function pointed to by *DataBlockReader with the parameter
                  BlockName.
                  When it finds data it calls the function pointed to by DataReader
                  with the parameter name and its value.
                  Returns 0 on success, -1 if file cannot be read, minus higher
                  number on other errors...

  Inputs  :
            int (*read_data_block) (char* block_name)
            int (*read_data) (char* par_name, char* par_value)
            char* file_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int read_init_file (int (*read_data_block) (char *block_name),
        int (*read_data) (char *par_name, char *par_value), char *file_name) {
    FILE *ini_file;     /* pointer to config file */
    int  i = 1;         /* to count lines */
    char *temp;         /* to read current line */
    char *temp_cut;     /* temp without its last letter */
    char *value;        /* value of current parameter */
    int  length;        /* to store string_lenth */
    int  j;             /* go through string */
    int  equ_sign = 0;  /* position of =-sign in line */
    int  result;        /* for storing returned values */

    ini_file = fopen(file_name, "rt");
    if (ini_file == NULL) {
        printf ("Error: File %s cannot be opened for reading.\n", file_name);
        return -4000;
    }

    /* reserve some space to read lines */
    temp = (char*) malloc (sizeof (char) * 512);
    if (temp == NULL) return -4001;
    temp_cut = (char*) malloc (sizeof (char) * 512);
    if (temp_cut==NULL)return -4002;
    value = (char*) malloc (sizeof (char) * 512);
    if (value == NULL) return -4003;

    while((i < 10000) && (!feof (ini_file))) {  /* read 10000 lines max or until file ends */
        if (fgets (temp, 511, ini_file) != NULL) {
            /* now the line has been read */
            /* check if first character is the comment sign # or lines is empty */
            /* then do not process the line */
            if ((temp[0] != '#') && (temp[0] != '\n')) {  /* line should be processed */
                length = (int) strlen (temp);

                if (length > 1) {
                    if (temp[0] == '[') {  /* check if first sign is [. then its a new block */
                        /* data Block handling */
                        /* extract the data block name */
                        j = 1;
                        while ((j < length) && (temp[j] != ']')) {
                            value[j-1] = temp[j];
                            j++;
                        }
                        value[j-1] = '\0';
                        result = read_data_block (value);
                        if (result != 0) return result;
                    }
                    else {  /* assume its a data line */
                        /* first, search for the "="-sign */
                        j = 0;
                        while (j < length) {
                           if (temp[j] == '=') {
                                equ_sign = j;
                                j = length;
                           }
                            j++;
                        }

                        /* separate the value from the current line assignment */
                        for (j=equ_sign+1; j<length; j++) {
                            value[j-equ_sign-1] = temp[j];

                            /* this should prevent some problems when reading
                               windows-generated config files in unix-like systems */
                            if (temp[j] == '\015') value[j-equ_sign-1] = '\0';
                        }
                        value[j-equ_sign-2] = '\0';

                        /* cut line before = sign */
                        temp[equ_sign] = '\0';
                        temp_cut[0] = '\0';

                        /* copy temp to temp_cut, but without last character */
                        for (j=0; j<(equ_sign-1); j++) temp_cut[j] = temp[j];
                        if (equ_sign > 0) temp_cut[equ_sign-1] = '\0';

                        /* let the designated function handle the data read */
                        result = read_data (temp, value);
                        if (result != 0) return result;
                    }
                }
            }
        }
    }

    fclose (ini_file);

    free (temp);
    free (temp_cut);
    free (value);

    return 0;
}

/*=============================================================================
  Function Name : write_int_array_to_file
  Description   : Writes the designated array of Count elements into a file.
                  If the file is just one column of values then FileType should
                  be set to 1. If the file contains 4 columns (x, y, z, value)
                  then set it to 0.

  Inputs  :
            char *file_name
            int  *source_array
            int  count
            int  file_type
  Outputs : no.

  Notes :
      The caller needs to make sure, that the array is large enough.

      Function call:
        target - int get_target_XYZ (int index, int* x, int* y, int* z);
=============================================================================*/
int write_int_array_to_file (char *file_name, int *source_array, int count, int file_type) {
    FILE *fp;
    int  i;
    int  x, y, z;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -4005;

    i = 0;
    switch (file_type) {
    case 1 :  /* just one column */
        if (normalize_output == 1) {
            for (i=0; i<count; i++)
                fprintf (fp, "%g\n", (double) source_array[i] * unit_conversion_factor);
        }
        else {
            for (i=0; i<count; i++) fprintf (fp, "%i\n", source_array[i]);
        }
        break;
    default :  /* four column file, x,y,z and value. */
       if (normalize_output == 1) {
            for (i=0; i<count; i++) {
                get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
                fprintf (fp, "%i\t%i\t%i\t%g\n", x, y, z, source_array[i] * unit_conversion_factor);
            }
        }
        else {
            for (i=0; i<count; i++) {
                get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
                fprintf (fp, "%i\t%i\t%i\t%i\n", x, y, z, source_array[i]);
            }
        }
        break;
    }

    fclose (fp);

    return 0;
}

/*=============================================================================
  Function Name : write_int_array_to_file
  Description   : Writes the designated array of Count elements into a file.
                  If the file is just one column of values then FileType should
                  be set to 1. If the file contains 4 columns (x, y, z, value)
                  then set it to 0.

  Inputs  :
            char* file_name
            int* source_array
            int count
            int file_type
  Outputs : no.

  Notes :
      The caller needs to make sure, that the array is large enough.

      Function call:
        target - int get_target_XYZ (int index, int* x, int* y, int* z);
=============================================================================*/
int write_double_array_to_file (char *file_name, double *source_array, int count, int file_type) {
    FILE *fp;
    int  i;
    int  x, y, z;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -4006;

    i = 0;
    switch (file_type) {
    case 1 :  /* just one column */
        for (i=0; i<count; i++) fprintf (fp, "%lg\n", source_array[i]);
        break;
    default :  /* four column file, x,y,z and value. */
        for (i=0; i<count; i++) {
            get_target_XYZ (i, &x, &y, &z);  /* convert linear index to 3 coords */
            fprintf (fp, "%i\t%i\t%i\t%lg\n", x, y, z, source_array[i] * unit_conversion_factor);
        }
        break;
    }

    fclose (fp);

    return 0;
}

/*=============================================================================
  Function Name : read_config_file_data_block
  Description   : Reads a data block from the configuration file.

  Inputs  :
            char* block_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int read_config_file_data_block (char *block_name) {
    /* ... we will ignore data blocks of the configuration file for now */
    return 0;
}

/*=============================================================================
  Function Name : read_config_file_data
  Description   : Reads data from the configuration file.
                  Needs to be called from the ini file reader while the
                  configuration input file is read.
                  Compare the parameter name to known parameters and then
                  read in corresponding value.

  Inputs  :
            char* par_name
            char* par_value
  Outputs : no.

  Notes : no.
=============================================================================*/
int read_config_file_data (char *par_name, char *par_value) {
    int itemp = 0;
    float ion_ax, ion_ay, ion_az;

    /* [IonBeam], ion beam parameters */
    if (strcmp (par_name, "max_no_ions") == 0){  /* read maximum number of ions to simulate */
        sscanf (par_value, "%i", &max_no_ions);
        if (print_level >= 1) printf ("Max number of ions:\t\t %i\n", max_no_ions);
    }
    if (strcmp (par_name, "ion_Z") == 0) {  /* read proton number of ion */
        sscanf (par_value, "%i", &ion_Z);
        if (print_level >= 1) printf ("ion beam Z:\t\t\t %i\n", ion_Z);
    }
    if (strcmp (par_name, "ion_M") == 0) {  /* read mass of ion */
        sscanf (par_value, "%f", &ion_M);
        if (print_level >= 1) printf ("ion beam M:\t\t\t %g amu\n", ion_M);
    }
    if (strcmp (par_name, "ion_E0") == 0) {  /* read mass of ion */
        sscanf (par_value, "%f", &ion_initial_energy);
        if (print_level >= 1) printf ("ion beam E_0:\t\t\t %g eV\n", ion_initial_energy);
    }
    if (strcmp (par_name, "ion_vx") == 0){  /* velocity vector of unit length */
        sscanf (par_value, "%f", &ion_vx);
        if (print_level >= 1) printf ("ion beam vx:\t\t\t %g\n", ion_vx);
    }
    if (strcmp (par_name, "ion_vy") == 0) {  /* velocity vector of unit length */
        sscanf (par_value, "%f", &ion_vy);
        if (print_level >= 1) printf ("ion beam vy:\t\t\t %g\n", ion_vy);
    }
    if (strcmp (par_name, "ion_vz") == 0) {  /* velocity vector of unit length */
        sscanf (par_value, "%f", &ion_vz);
        if (print_level >= 1) printf ("ion beam vz:\t\t\t %g\n", ion_vz);
    }
    if (strcmp (par_name, "ion_ax") == 0){  /* velocity vector of unit length */
        sscanf (par_value, "%f", &ion_ax);
        ion_vx = cos (ion_ax * D_TO_R);
        if (print_level >= 1) printf ("ion beam vx:\t\t\t %g\n", ion_vx);
    }
    if (strcmp (par_name, "ion_ay") == 0) {  /* velocity vector of unit length */
        sscanf (par_value, "%f", &ion_ay);
        ion_vy = cos (ion_ay * D_TO_R);
        if (print_level >= 1) printf ("ion beam vy:\t\t\t %g\n", ion_vy);
    }
    if (strcmp (par_name, "ion_az") == 0) {  /* velocity vector of unit length */
        sscanf (par_value, "%f", &ion_az);
        ion_vz = cos (ion_az * D_TO_R);
        if (print_level >= 1) printf ("ion beam vz:\t\t\t %g\n", ion_vz);
    }
    if (strcmp (par_name, "ion_distribution") == 0) {  /* ... */
        sscanf (par_value, "%i", &ion_distribution);
        if (print_level >= 1) printf ("Ion entry distribution model:\t %i.\n", ion_distribution);
    }
    if (strcmp (par_name, "enter_x") == 0) {  /* point of entry */
        sscanf (par_value,"%f", &enter_x);
        if (print_level >= 1) printf ("Entry point x:\t\t\t %g nm\n", enter_x);
    }
    if (strcmp (par_name, "enter_y") == 0) {  /* point of entry */
        sscanf (par_value,"%f", &enter_y);
        if (print_level >= 1) printf ("Entry point y:\t\t\t %g nm\n", enter_y);
    }
    if (strcmp (par_name, "enter_z") == 0) {  /* point of entry */
        sscanf (par_value, "%f", &enter_z);
        if (print_level >= 1) printf ("Entry point z:\t\t\t %g nm\n", enter_z);
    }
    if (strcmp (par_name, "beam_spread") == 0) {
        sscanf (par_value, "%f", &beam_spread);
        if (print_level>=1) printf ("Beam spread:\t\t\t %g nm\n", beam_spread);
    }

    /* for compatibility reasons: */
    if (strcmp (par_name, "ion_angle_y") == 0) {  /* read orientation of ion beam */
        printf ("ERROR! ion_angle_y option has been disabled!\n");
        return -4007;
    }
    /* for compatibility reasons: */
    if (strcmp (par_name, "ion_angle_z") == 0) {  /* read orientation of ion beam */
        printf ("ERROR! ion_angle_z option has been disabled!\n");
        return -4008;
    }

    /* [Simulation] */
    if (strcmp (par_name, "OutputFileBaseName") == 0) {  /* ... */
        OutputFileBaseName = malloc (sizeof (char) * (strlen (par_value) + 1));
        if (OutputFileBaseName == NULL) return -4009;
        strcpy (OutputFileBaseName, par_value);
    }
    if (strcmp (par_name, "output_format") == 0) {  /* ... */
        sscanf (par_value, "%i", &output_format);
        if (print_level >= 1) printf ("Output format :\t\t %i \n", output_format);
    }
    if (strcmp (par_name, "normalize_output") == 0) {  /* ... */
        sscanf (par_value, "%i", &normalize_output);
        if (print_level >= 1) printf ("Normalize Output :\t\t %i \n", normalize_output);
    }
    if (strcmp (par_name, "display_interval") == 0){  /* ... */
        sscanf (par_value, "%i", &display_interval);
        if (print_level >= 1) printf ("Show ion number every \t\t %i ions.\n", display_interval);
    }
    if (strcmp (par_name, "storage_interval") == 0) {  /* ... */
        sscanf (par_value, "%i", &storage_interval);
        if (print_level >= 1) printf ("Store results every \t\t %i ions.\n", storage_interval);
    }
    if (strcmp (par_name, "status_update_interval") == 0){  /* ... */
        sscanf (par_value, "%i", &itemp);
        if (itemp>0) status_update_interval = itemp;
        if (print_level >= 1) printf ("Write status file every \t\t %i ions.\n", status_update_interval);
    }
    if (strcmp (par_name, "store_transmitted_ions") == 0) {  /* store transmitted ions? */
        sscanf (par_value, "%i", &store_transmitted_ions);
       if (print_level >= 1) printf ("Storing transmitted ions?\t %i\n", store_transmitted_ions);
    }
    if (strcmp (par_name, "store_exiting_recoils") == 0) {  /* store exiting recoils? */
        sscanf (par_value, "%i", &store_exiting_recoils);
        if (print_level >= 1) printf ("Storing exiting recoils?\t %i\n", store_exiting_recoils);
    }
    if (strcmp (par_name, "store_exiting_limit") == 0) {  /* store transmitted recoils until this number */
        sscanf (par_value, "%i", &store_exiting_limit);
        if (print_level >= 1) printf ("Max. exiting recoils:\t\t%i\n", store_exiting_limit);
    }
    /* if energy deposition should be stored, this must be 1 */
    if (strcmp (par_name, "store_energy_deposit") == 0) {
        sscanf (par_value, "%i", &store_energy_deposit);
        if (print_level >= 1) printf ("Detailed energy deposition:\t %i\n", store_energy_deposit);
    }
    if (strcmp (par_name, "store_ion_paths") == 0) {  /* store ion paths? */
        sscanf (par_value, "%i", &store_ion_paths);
        if (print_level >= 1) printf ("Storing ions paths?\t\t %i\n", store_ion_paths);
    }
    if (strcmp (par_name, "store_recoil_cascades") == 0) {  /* store recoil cascades? */
        sscanf (par_value, "%i", &store_recoil_cascades);
        if (print_level >= 1) printf ("Storing exact recoils cascades?\t %i\n", store_recoil_cascades);
    }
    if (strcmp (par_name, "store_path_limit") == 0) {  /* ... */
        sscanf (par_value, "%i", &store_path_limit);
        if (print_level >= 1) printf ("Store paths or cascades only until ion %i.\n", store_path_limit);
    }
    /* how detailed the simulation should be */
    if (strcmp (par_name, "simulation_type") == 0) {
        sscanf (par_value, "%i", &simulation_type);
        if (print_level >= 1) printf ("Simulation type:\t\t %i\n", simulation_type);
    }
    if (strcmp (par_name, "transport_type") == 0) {  /* ... */
        sscanf (par_value, "%i", &transport_type);
        if (print_level >= 1) printf ("Transport type:\t\t\t %i \n", transport_type);
    }
    if (strcmp (par_name, "flight_length_type") == 0) {  /* ... */
        sscanf (par_value, "%i", &flight_length_type);
        if (print_level >= 1) printf ("Flight length distribution model:%i\n", flight_length_type);
    }
    if (strcmp (par_name, "flight_length_constant") == 0) {  /* ... */
        sscanf (par_value, "%f", &flight_length_const);
        if (print_level >= 1) printf ("Flight length:\t\t\t %g nm\n", flight_length_const);
    }
    if (strcmp (par_name, "scattering_calculation") == 0) {  /* ... */
        sscanf (par_value, "%i", &scattering_calculation);
        if (print_level >= 1) printf ("Scattering calculation:\t\t %i \n", scattering_calculation);
    }
    if (strcmp (par_name, "tracing_recoil_or_not") == 0) {  /* KP model, Aug. 5, 2014 */
        sscanf (par_value, "%i", &tracing_recoil_or_not);
        if (print_level >= 1) printf ("Tracing recoil or not:\t\t %i \n", tracing_recoil_or_not);
    }
    if (strcmp (par_name, "multiple_collisions") == 0){  /* ... */
        sscanf (par_value, "%i", &max_annular_coll_volumes);
        if (print_level >= 1) printf ("Multiple collisions :\t\t %i \n", max_annular_coll_volumes);
    }
    /* how detailed the sputtering should be evaluated */
    if (strcmp (par_name, "detailed_sputtering") == 0) {
        sscanf (par_value, "%i", &detailed_sputtering);
        if (print_level >= 1) printf ("Detailed sputtering:\t\t %i\n", detailed_sputtering);
    }
    if (strcmp (par_name, "single_ion_sputter_yields") == 0) {  /* ... */
        sscanf (par_value, "%i", &single_ion_sputter_yields);
        if (print_level >= 1) printf ("single_ion_sputter_yields:\t %i \n", single_ion_sputter_yields);
    }
    if (strcmp (par_name, "do_not_store_damage") == 0) {  /* ... */
        sscanf (par_value, "%i", &do_not_store_damage);
        if (print_level >= 1) printf ("Do not store damage:\t\t %i \n", do_not_store_damage);
    }
    /* minimum energy for moving projectiles */
    if (strcmp (par_name, "min_energy") == 0) {
        sscanf (par_value, "%f", &min_energy);
        if (print_level >= 1) printf ("Minimum projectile energy:\t %g eV\n", min_energy);
    }
    if (strcmp (par_name, "seed1") == 0) {  /* ... */
        sscanf (par_value, "%i", &seed1);
        if (print_level >= 1) printf ("Random seed 1 set to:\t\t %i.\n", seed1);
    }
    if (strcmp (par_name, "seed2") == 0){  /* ... */
        sscanf (par_value, "%i", &seed2);
        if (print_level >= 1) printf ("Random seed 2 set to:\t\t %i.\n", seed2);
    }

    /* [Target] */
    if (strcmp (par_name, "geometry_type") == 0) {
        sscanf (par_value, "%i", &geometry_type);
        if (print_level >= 1) printf ("geometry_type:\t\t\t %i \n", geometry_type);
    }
    if (strcmp (par_name, "no_substrate") == 0) {
        sscanf (par_value, "%i", &no_substrate);
        if (print_level >= 1) printf ("no_substrate:\t\t\t %i \n", no_substrate);
    }
    if (strcmp (par_name, "gen_shape_or_not") == 0) {
        sscanf (par_value, "%i", &gen_shape_or_not);
        if (print_level >= 1) printf ("gen_shape_or_not:\t\t\t %i \n", gen_shape_or_not);
    }
    if (strcmp (par_name, "straggling_model") == 0) {  /* ... */
        sscanf (par_value, "%i", &straggling_model);
        if (print_level >= 1) printf ("Straggling model is:\t\t %i\n", straggling_model);
    }
    if (strcmp (par_name, "MaterialsFileName") == 0) {  /* ... */
        if (single_input_file != 1) strcpy (MaterialsFileName, par_value);
    }
    if (strcmp (par_name, "TargetstructureFileName") == 0) {  /* ... */
        if (single_input_file != 1) strcpy (TargetStructureFileName, par_value);
    }

    return 0;
}

/*=============================================================================
  Function Name : store_results_iradina
  Description   : Store the results of the simulation (arrays with distribution
                  of implanted ions, defects etc.)
                  This is the non-dynamic version of the function for
                  material-based output.

  Inputs  :
            char* base_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int store_results_iradina (char *base_name) {
    int  base_name_length;
    char *str_temp;
    char *str_temp2;
    int i,j;
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
    strcat (str_temp, ".ions.total");
    if (print_level >= 2) printf ("Storing implanted ions to:      %s\n", str_temp);
    write_int_array_to_file (str_temp, target_implanted_ions, cell_count, target_composition_file_type);
    str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */

    if (do_not_store_damage == 0) {
        /* Part of ions that replaced identical target atoms */
        strcat (str_temp, ".ions.replacements");
        if (print_level >= 2) printf ("Storing implanted replacing ions to:   %s\n", str_temp);
        write_int_array_to_file (str_temp, target_replacing_ions, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    }

    sum_up_material_arrays ();  /* sum up the individual element-dependent arrays */

    if (do_not_store_damage == 0) {
        /* sum of vacancies: */
        strcat (str_temp, ".vac.sum");
        if (print_level>=2) printf ("Storing total vacancies to: %s\n", str_temp);
        write_int_array_to_file (str_temp, target_total_vacancies, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';

        /* sum of displacements: */
        strcat (str_temp, ".disp.sum");
        if (print_level >= 2) printf ("Storing total displacements to: %s\n", str_temp);
        write_int_array_to_file (str_temp, target_total_displacements, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';

        /* sum of recoiled interstitials: */
        strcat (str_temp, ".int.sum");
        if (print_level >= 2) printf ("Storing sum of recoiled interstitials:   %s\n", str_temp);
        write_int_array_to_file (str_temp, target_total_interstitials, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';

        /* sum of recoil replacements: */
        strcat (str_temp, ".repl.sum");
        if (print_level >= 2) printf ("Storing sum of recoiled replacements:   %s\n", str_temp);
        write_int_array_to_file (str_temp, target_total_replacements, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';
    }

    /* deposited energy */
    if (store_energy_deposit == 1) {
        strcat (str_temp, ".energy.phonons");
        if (print_level >= 2) printf ("Storing energy deposited to phonons:  %s\n", str_temp);
        write_double_array_to_file (str_temp, target_energy_phonons, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';
        strcat (str_temp, ".energy.electronic");
        if (print_level >= 2) printf ("Storing energy deposited to electrons:%s\n", str_temp);
        write_double_array_to_file (str_temp, target_energy_electrons, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';
    }

    /* store single ion sputter yields */
    if (single_ion_sputter_yields == 1) {
        strcat (str_temp, ".single_ion_sputter_yields");
        if (print_level >= 2) printf ("Storing single ion sputter yields to:   %s\n", str_temp);
        fp = fopen (str_temp, "w");
        if (fp == NULL) {
            printf("Error: cannot store histogram.\n");
            return -4010;
        }
        for (i=0; i<=MAX_SPUTTERED; i++) fprintf (fp, "%i\t%i\n", i, sputter_yield_histogram[i]);
        fclose (fp);
        str_temp[base_name_length] = '\0';
    }

    /* sum of atoms leaving the target */
    if (detailed_sputtering == 1) {
        /* particles leaving the sample (sputtered or implanted deeper): */
        strcat (str_temp, ".leaving_directions.sum");
        if (print_level >= 2) printf (" Storing sum of leaving atoms to:       %s\n", str_temp);
        sprintf (str_temp2, "(+x,+y,+z)\t%i\n(+x,-y,+z)\t%i\n(+x,-y,-z)\t%i\n(+x,+y,-z)\t%i\n(-x,+y,+z)\t%i\n(-x,-y,+z)\t%i\n(-x,-y,-z)\t%i\n(-x,+y,-z)\t%i\n",
                 total_sputter_counter[0], total_sputter_counter[1], total_sputter_counter[2], total_sputter_counter[3],
                 total_sputter_counter[4], total_sputter_counter[5], total_sputter_counter[6], total_sputter_counter[7]);
        write_string_to_file (str_temp, str_temp2);
        str_temp [base_name_length] = '\0';

        strcat (str_temp, ".leaving.sum");
        if (print_level >= 2) printf ("Storing sum of leaving atoms (cells):   %s\n", str_temp);
        write_int_array_to_file (str_temp, target_total_sputtered, cell_count, target_composition_file_type);
        str_temp[base_name_length] = '\0';

        /* ions leaving the target */
        strcat (str_temp, ".leaving_directions.ions");
        if (print_level >= 2) printf (" Storing sum of leaving ions to:        %s\n", str_temp);
        sprintf (str_temp2, "(+x,+y,+z)\t%i\n(+x,-y,+z)\t%i\n(+x,-y,-z)\t%i\n(+x,+y,-z)\t%i\n(-x,+y,+z)\t%i\n(-x,-y,+z)\t%i\n(-x,-y,-z)\t%i\n(-x,+y,-z)\t%i\n",
                 leaving_ions[0], leaving_ions[1], leaving_ions[2], leaving_ions[3],
                 leaving_ions[4], leaving_ions[5], leaving_ions[6], leaving_ions[7]);
        write_string_to_file (str_temp, str_temp2);
        str_temp[base_name_length] = '\0';
    }

    /* sum of ints and vacs and so on for each element in each material */
    for (i=0; i<number_of_materials; i++) {
        for (j=0; j<list_of_materials[i].element_count; j++) {  /* Go through elements, store arrays for each */
            if (print_level >= 2) printf ("Storing for material %i, elem %i:\n", i, j);

            if (do_not_store_damage == 0) {
                sprintf (str_temp+base_name_length, ".int.z%i.m%.3f.mat%i.elem%i",
                         list_of_materials[i].elements_Z[j],
                         list_of_materials[i].elements_M[j], i, j);
                if (print_level >= 2) printf (" Recoil interstitials to %s\n", str_temp);
                write_int_array_to_file(str_temp, list_of_materials[i].target_implanted_recoils_int[j],
                                        cell_count, target_composition_file_type);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, ".repl.z%i.m%.3f.mat%i.elem%i",
                         list_of_materials[i].elements_Z[j],
                         list_of_materials[i].elements_M[j], i, j);
                if (print_level >= 2) printf (" Recoil replacements to %s\n", str_temp);
                write_int_array_to_file (str_temp, list_of_materials[i].target_implanted_recoils_repl[j],
                                         cell_count, target_composition_file_type);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, ".vac.z%i.m%.3f.mat%i.elem%i",
                         list_of_materials[i].elements_Z[j],
                         list_of_materials[i].elements_M[j], i, j);
                if (print_level >= 2) printf (" Vacancies to     %s\n", str_temp);
                write_int_array_to_file (str_temp, list_of_materials[i].target_elemental_vacancies[j],
                                         cell_count, target_composition_file_type);
                str_temp[base_name_length] = '\0';

                sprintf (str_temp+base_name_length, ".disp.z%i.m%.3f.mat%i.elem%i",
                         list_of_materials[i].elements_Z[j],
                         list_of_materials[i].elements_M[j], i, j);
                if (print_level >= 2) printf (" Displacements to   %s\n", str_temp);
                write_int_array_to_file (str_temp, list_of_materials[i].target_elemental_disp[j],
                                         cell_count, target_composition_file_type);
                str_temp[base_name_length] = '\0';
            }

            if (detailed_sputtering == 1) {
                /* particles leaving the sample (sputtered or implanted deeper): */
                sprintf (str_temp+base_name_length, ".leaving_directions.z%i.m%.3f.mat%i.elem%i",
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

                sprintf (str_temp+base_name_length, ".leaving.z%i.m%.3f.mat%i.elem%i",
                         list_of_materials[i].elements_Z[j],
                         list_of_materials[i].elements_M[j], i, j);
                if (print_level >= 2) printf (" Leaving from cell  %s\n", str_temp);
                write_int_array_to_file (str_temp, list_of_materials[i].target_sputtered_atoms[j],
                                         cell_count, target_composition_file_type);
                str_temp[base_name_length] = '\0';
            }
        }
    }

    if (store_transmitted_ions == 1) {
        strcat (str_temp, ".transmitted.ions");
        if (print_level >= 2) printf ("Storing transmitted ions to: %s\n", str_temp);
        store_transmission_array (str_temp, transmit_list, transmission_pointer);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    }
    if (store_exiting_recoils == 1) {
        if (print_level >= 2) printf ("Storing transmitted recoils...\n");
        for (i=0; i<number_of_materials; i++) {  /* go through mats */
            for (j=0; j<list_of_materials[i].element_count; j++) { /* go through elements, store arrays for each */
                sprintf (str_temp+base_name_length, ".leaving_recoils.z%i.m%.3f.mat%i.elem%i",
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
        strcat (str_temp, ".depth_dist_functions.dat");
        store_depth_dist_array (str_temp);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */

        do_radial_dist_statistics ();
        strcat (str_temp, ".radial_dist_functions.dat");
        store_radial_dist_array (str_temp);
        str_temp[base_name_length] = '\0';  /* reset strTemporary filename for next storing operation */
    //}

    free (str_temp);

    return 0;
}

/*=============================================================================
  Function Name : sum_up_material_arrays
  Description   : Interstitials and so on are stored for each element from each
                  material separately, but may also be interesting in sum.
                  So this function does all the summing up.

  Inputs  : no.
  Outputs : no.

  Notes :
      Function call:
        utils - int fill_zero (int* array, int count);
        void add_int_array (int* dest, int* source, int count).
=============================================================================*/
void sum_up_material_arrays (void) {
    int i, j, k;

    /* reset sum arrays: */
    fill_int_zero (target_total_vacancies    , cell_count);
    fill_int_zero (target_total_displacements, cell_count);
    fill_int_zero (target_total_interstitials, cell_count);
    fill_int_zero (target_total_replacements , cell_count);

    if (detailed_sputtering == 1) {
        /* arrays of leaving atoms */
        for (k=0; k<8; k++) total_sputter_counter[k] = 0;  /* reset sum */

        fill_int_zero (target_total_sputtered, cell_count);
    }

    for (i=0; i<number_of_materials; i++) {  /* loop through materials */
        for (j=0; j<list_of_materials[i].element_count; j++) {  /* loop through elements */
            add_int_array (target_total_vacancies,
                (list_of_materials[i].target_elemental_vacancies)[j]   , cell_count);
            add_int_array (target_total_replacements,
                (list_of_materials[i].target_implanted_recoils_repl)[j], cell_count);
            add_int_array (target_total_displacements,
                (list_of_materials[i].target_elemental_disp)[j]        , cell_count);
            add_int_array (target_total_interstitials,
                (list_of_materials[i].target_implanted_recoils_int)[j] , cell_count);
            if (detailed_sputtering == 1) {
                add_int_array (target_total_sputtered,
                    (list_of_materials[i].target_sputtered_atoms)[j]   , cell_count);

                for (k=0; k<8; k++)  /* loop through all 8 directions/quadrants */
                    total_sputter_counter[k] += list_of_materials[i].sputter_counter[(8*j)+k];
            }
        }
    }

    return;
}

/*=============================================================================
  Function Name : do_depth_dist_statistics
  Description   : Interstitials and so on are stored for each element from each
                  material separately, but may also be interesting in doing the
                  depth distribution summing up.
                  So this function gives the depth distributions.

  Inputs  : no.
  Outputs : no.

  Notes :
      Function call:
        utils - int fill_zero (int* array, int count);
        void add_int_array (int* dest, int* source, int count).
=============================================================================*/
void do_depth_dist_statistics (void) {
    int i, j, k;
    int x, y, z;

    /* reset depth dist arrays */
    fill_int_zero (target_depth_implanted_ions, cell_count_z);
    fill_int_zero (target_depth_replacing_ions, cell_count_z);
    fill_double_zero (target_depth_energy_electrons, cell_count_z);
    fill_double_zero (target_depth_energy_phonons, cell_count_z);

    fill_int_zero (target_depth_total_interstitials, cell_count_z);
    fill_int_zero (target_depth_total_replacements, cell_count_z);
    fill_int_zero (target_depth_total_vacancies, cell_count_z);
    fill_int_zero (target_depth_total_displacements, cell_count_z);

    for (i=0; i<number_of_materials; i++) {
        for (j=0; j<list_of_materials[i].element_count; j++) {
            fill_int_zero (list_of_materials[i].target_depth_implanted_recoils_int[j], cell_count_z);
            fill_int_zero (list_of_materials[i].target_depth_implanted_recoils_repl[j], cell_count_z);
            fill_int_zero (list_of_materials[i].target_depth_elemental_vacancies[j], cell_count_z);
            fill_int_zero (list_of_materials[i].target_depth_elemental_disp[j], cell_count_z);
        }
    }

    for (i=0; i<cell_count; i++) {
        /* convert linear index to 3 coords */
        get_target_XYZ (i, &x, &y, &z);

        target_depth_implanted_ions[z] += target_implanted_ions[i];
        target_depth_replacing_ions[z] += target_replacing_ions[i];
        target_depth_energy_electrons[z] += target_energy_electrons[i];
        target_depth_energy_phonons[z] += target_energy_phonons[i];

        target_depth_total_interstitials[z] += target_total_interstitials[i];
        target_depth_total_replacements[z] += target_total_replacements[i];
        target_depth_total_vacancies[z] += target_total_vacancies[i];
        target_depth_total_displacements[z] += target_total_displacements[i];

        for (j=0; j<number_of_materials; j++) {
            for (k=0; k<list_of_materials[j].element_count; k++) {
                list_of_materials[j].target_depth_implanted_recoils_int[k][z] +=
                    list_of_materials[j].target_implanted_recoils_int[k][i];
                list_of_materials[j].target_depth_implanted_recoils_repl[k][z] +=
                    list_of_materials[j].target_implanted_recoils_repl[k][i];
                list_of_materials[j].target_depth_elemental_vacancies[k][z] +=
                    list_of_materials[j].target_elemental_vacancies[k][i];
                list_of_materials[j].target_depth_elemental_disp[k][z] +=
                    list_of_materials[j].target_elemental_disp[k][i];
            }
        }
    }

    return;
}

/*=============================================================================
  Function Name : do_radial_dist_statistics
  Description   : Interstitials and so on are stored for each element from each
                  material separately, but may also be interesting in doing the
                  radial distribution summing up.
                  So this function gives the radial distributions.

  Inputs  : no.
  Outputs : no.

  Notes :
      Function call:
        utils - int fill_zero (int* array, int count);
        void add_int_array (int* dest, int* source, int count).
=============================================================================*/
void do_radial_dist_statistics (void) {
    int i, j, k;
    int x, y, z, ir;

    /* reset depth dist arrays */
    fill_int_zero (target_radial_implanted_ions, cell_max_xy);
    fill_int_zero (target_radial_replacing_ions, cell_max_xy);
    fill_double_zero (target_radial_energy_electrons, cell_max_xy);
    fill_double_zero (target_radial_energy_phonons, cell_max_xy);

    fill_int_zero (target_radial_total_interstitials, cell_max_xy);
    fill_int_zero (target_radial_total_replacements, cell_max_xy);
    fill_int_zero (target_radial_total_vacancies, cell_max_xy);
    fill_int_zero (target_radial_total_displacements, cell_max_xy);

    for (i=0; i<number_of_materials; i++) {
        for (j=0; j<list_of_materials[i].element_count; j++) {
            fill_int_zero (list_of_materials[i].target_radial_implanted_recoils_int[j], cell_max_xy);
            fill_int_zero (list_of_materials[i].target_radial_implanted_recoils_repl[j], cell_max_xy);
            fill_int_zero (list_of_materials[i].target_radial_elemental_vacancies[j], cell_max_xy);
            fill_int_zero (list_of_materials[i].target_radial_elemental_disp[j], cell_max_xy);
        }
    }

    for (i=0; i<cell_count; i++) {
        /* convert linear index to 3 coords */
        get_target_XYZ (i, &x, &y, &z);

        ir = (int) (sqrt ((float) (x * x + y * y)));

        target_radial_implanted_ions[ir] += target_implanted_ions[i];
        target_radial_replacing_ions[ir] += target_replacing_ions[i];
        target_radial_energy_electrons[ir] += target_energy_electrons[i];
        target_radial_energy_phonons[ir] += target_energy_phonons[i];

        target_radial_total_interstitials[ir] += target_total_interstitials[i];
        target_radial_total_replacements[ir] += target_total_replacements[i];
        target_radial_total_vacancies[ir] += target_total_vacancies[i];
        target_radial_total_displacements[ir] += target_total_displacements[i];

        for (j=0; j<number_of_materials; j++) {
            for (k=0; k<list_of_materials[j].element_count; k++) {
                list_of_materials[j].target_radial_implanted_recoils_int[k][ir] +=
                    list_of_materials[j].target_implanted_recoils_int[k][i];
                list_of_materials[j].target_radial_implanted_recoils_repl[k][ir] +=
                    list_of_materials[j].target_implanted_recoils_repl[k][i];
                list_of_materials[j].target_radial_elemental_vacancies[k][ir] +=
                    list_of_materials[j].target_elemental_vacancies[k][i];
                list_of_materials[j].target_radial_elemental_disp[k][ir] +=
                    list_of_materials[j].target_elemental_disp[k][i];
            }
        }
    }

    return;
}

/*=============================================================================
  Function Name : store_depth_dist_array
  Description   : Store array of depth distribution functions.

  Inputs  :
            char* file_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int store_depth_dist_array (char *file_name) {
    FILE *fp;
    int i;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -4011;

    for (i=0; i<cell_count_z; i++)
        fprintf (fp, "%i\t%i\t%i\t%g\t%g\t%i\t%i\t%i\t%i\n", i,
                 target_depth_implanted_ions[i],      target_depth_replacing_ions[i],
                 target_depth_energy_electrons[i],    target_depth_energy_phonons[i],
                 target_depth_total_interstitials[i], target_depth_total_replacements[i],
                 target_depth_total_vacancies[i],     target_depth_total_displacements[i]);

    //TODO: list_of_materials[i].target_depth_*[j]

    return 0;
}

/*=============================================================================
  Function Name : store_radial_dist_array
  Description   : Store array of radial distribution functions.

  Inputs  :
            char* file_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int store_radial_dist_array (char *file_name) {
    FILE *fp;
    int i;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -4011;

    for (i=0; i<cell_max_xy; i++)
        fprintf (fp, "%i\t%i\t%i\t%g\t%g\t%i\t%i\t%i\t%i\n", i,
                 target_radial_implanted_ions[i],      target_radial_replacing_ions[i],
                 target_radial_energy_electrons[i],    target_radial_energy_phonons[i],
                 target_radial_total_interstitials[i], target_radial_total_replacements[i],
                 target_radial_total_vacancies[i],     target_radial_total_displacements[i]);

    //TODO: list_of_materials[i].target_radial_*[j]


    return 0;
}

/*=============================================================================
  Function Name : store_transmission_array
  Description   : Store array of transmitted ions.

  Inputs  :
            char* file_name
            struct transmitted_ion* trans_array
            int tr_pointer
  Outputs : no.

  Notes : no.
=============================================================================*/
int store_transmission_array (char *file_name, struct transmitted_ion *trans_array,
                              int tr_pointer) {
    FILE *fp;
    int  i;

    fp = fopen (file_name, "wt");
    if (fp == NULL) return -4012;

    i = 0;
    for (i=0; i<tr_pointer; i++)
        fprintf (fp, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                 trans_array[i].x,
                 trans_array[i].y,
                 trans_array[i].z,
                 trans_array[i].vx,
                 trans_array[i].vy,
                 trans_array[i].vz,
                 trans_array[i].energy);

    fclose (fp);

    return 0;
}

/*=============================================================================
  Function Name : open_file_continuous
  Description   : Opens file base_name + extension, keeps it open.

  Inputs  :
            char* base_name
            char* extension
  Outputs : no.

  Notes : no.
=============================================================================*/
FILE *open_file_continuous (char *base_name, char *extension) {
    int  base_name_length;
    char *strTemp;

    base_name_length = strlen (base_name);
    strTemp = (char*) malloc (sizeof (char) * (base_name_length + 50));
    strncpy (strTemp, base_name, base_name_length+1);
    strcat  (strTemp, extension);

    return fopen (strTemp, "wt");
}

/*=============================================================================
  Function Name : read_float_block
  Description   : Reads a block of Count float values from file file_name
                  starting at Offset and puts these in Array.

  Inputs  :
            char* file_name
            int offset
  Outputs :
            int count
            float* array

  Notes :
      Might cause a seg fault or crash if array is too short!
=============================================================================*/
int read_float_block (char *file_name, int offset, int count, float *array) {
    int iscan;
    FILE *fp;

    fp = fopen (file_name, "rb");
    if (fp == NULL) return -4013;  /* cannot open file */
    fseek (fp, sizeof (float) * offset, 1);
    iscan = fread ((void*) array, sizeof (float), count, fp);  /* read data block */
    if (iscan < 0) return -4014;
    fclose (fp);

    return 0;
}

/*=============================================================================
  Function Name : write_string_to_file
  Description   : Does what you think it does.
                  Returns 0 on success.

  Inputs  :
            char* file_name
            char* str
  Outputs : no.

  Notes : no.
=============================================================================*/
int write_string_to_file (char *file_name, char *str) {
    FILE *fp;

    fp = fopen (file_name, "w");
    if (fp == NULL) {
        return -4015;
    }
    fprintf (fp, "%s", str);
    fclose (fp);

    return 0;
}

/*=============================================================================
  Function Name : dispaly_a_file
  Description   : Prints the contents of a text file to the std out.
                  Returns 0 on success.

  Inputs  : char* Filename.
  Outputs : no.

  Notes : no.
=============================================================================*/
int display_a_file (char *file_name) {
    FILE *fp;

    fp = fopen (file_name, "r");
    if (fp == NULL){
        return -4016;
    }
    while (!feof(fp)) {
        putchar (fgetc (fp));
    }

    return 0;
}

/*=============================================================================
  Function Name : split_single_input_file
  Description   : If a single input file is provided, split it to temp files
                  that can be read later.

  Inputs  : char* file_name - the name of the combined input file.
  Outputs : no.

  Notes : no.
=============================================================================*/
int split_single_input_file (char *file_name) {
    FILE *fp;       /* pointer to combined input file */
    FILE *fout[4];  /* pointer to the four single output files */
    int  i;
    int  file_num;   /* file, into which is currently written */
    int  sep_line;   /* if the current line separates files, this is 1 */
    char *line;     /* to read current line */
    char *temp_FN;  /* to store filename */

    /* make a copy of the config filename, so we can change it later */
    temp_FN = (char*) malloc (sizeof (char) * 1024);
    if (temp_FN == NULL) {printf ("Error. Insufficient memory!\n"); return -4017;}
    strncpy (temp_FN, ConfigFileName, 1023);

    /* name of the general input config file */
    strcpy (ConfigFileName, "temp_configfile.im3d");
    /* filename of the file that define the structure of the target */
    strcpy (TargetStructureFileName, "temp_structfile.im3d");
    /* target composition file */
    strcpy (TargetCompositionFileName, "temp_compfile.im3d");
    /* name of the file that defines the materials in the target */
    strcpy (MaterialsFileName, "temp_matfile.im3d");

    /* remove existing temp files: */
    remove (ConfigFileName);
    remove (TargetStructureFileName);
    remove (TargetCompositionFileName);
    remove (MaterialsFileName);
    single_input_file = 1;  /* mark that this is used indeed */

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=1) printf ("Splitting combined input file...\n");
    /* MPI=============================================== */
#else
    if (print_level >= 1) printf ("Splitting combined input file...\n");
#endif

    /* open input file */
    fp = fopen (temp_FN, "r");
    if (fp == NULL) {
        printf ("Error: Cannot open combined config file: %s\n", file_name);
        return -4018;
    }

    /* open output files */
    fout[0] = fopen (ConfigFileName, "w");
    if (fout[0] == NULL) {
        printf ("Error: Cannot open file for writing: %s\n", ConfigFileName);
        return -4019;
    }
    fout[1] = fopen (TargetStructureFileName, "w");
    if (fout[1] == NULL) {
        printf ("Error: Cannot open file for writing: %s\n", TargetStructureFileName);
        return -4020;
    }
    fout[2] = fopen (MaterialsFileName, "w");
    if (fout[2] == NULL) {
        printf ("Error: Cannot open file for writing: %s\n", MaterialsFileName);
        return -4021;
    }
    fout[3] = fopen (TargetCompositionFileName, "w");
    if (fout[3] == NULL) {
        printf ("Error: Cannot open file for writing: %s\n", TargetCompositionFileName);
        return -4022;
    }

    line = (char*) malloc (sizeof (char) * 512);  /* reserve space to read line */
    if (line == NULL) return -4023;  /* out of mem */

    file_num = 0;
    while (!feof(fp)) {  /* until file ends */
        sep_line = 0;
        if (fgets (line, 511, fp) != NULL) {  /* now the line has been read */
            /* check for file separators */
            if (strncmp (line, "#<<<BEGIN CONFIGFILE", 20) == 0) {  /* it's the config
                                                                    file that follows */
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node==ROOT && print_level>=2) printf ("Config file part found.\n");
                /* MPI=============================================== */
#else
                if (print_level >= 2) printf ("Config file part found.\n");
#endif
                file_num = 0;
                sep_line = 1;
            }
            if (strncmp (line, "#<<<BEGIN STRUCTUREFILE", 23) == 0){  /* it's the structures
                                                                        file that follows */
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node==ROOT && print_level>=2) printf ("Structure file part found.\n");
                /* MPI=============================================== */
#else
                if (print_level >= 2) printf ("Structure file part found.\n");
#endif
                file_num = 1;
                sep_line = 1;
            }
            if (strncmp (line, "#<<<BEGIN MATFILE", 17) == 0) {  /* it's the material file
                                                                    that follows */
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node==ROOT && print_level>=2) printf ("Material file part found.\n");
                /* MPI=============================================== */
#else
                if (print_level >= 2) printf ("Material file part found.\n");
#endif
                file_num = 2;
                sep_line = 1;
            }
            if (strncmp (line, "#<<<BEGIN COMPFILE",18) == 0) {  /* it's the composition
                                                                    file that follows */
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node==ROOT && print_level>=2) printf ("Composition file part found.\n");
                /* MPI=============================================== */
#else
                if (print_level >= 2) printf ("Composition file part found.\n");
#endif
                file_num = 3;
                sep_line = 1;
            }

            if (sep_line == 0) {  /* its a data line, put it into current output file */
                fputs (line, fout[file_num]);
            }
        }
    }

    /* close files */
    fclose (fp);
    for (i=0; i<4; i++) fclose (fout[i]);

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && print_level>=1) printf ("Splitting finished successfully.\n");
    /* MPI=============================================== */
#else
    if (print_level >= 1) printf ("Splitting finished successfully.\n");
#endif

    return 0;
}

/*=============================================================================
  Function Name : check_split_input_file
  Description   : Checks if the given config file is a single combined input
                  file (incl. structure, materials, and composition).
                  If so, the file is split into four separate temp files for
                  further "conventional" processing.

  Inputs  : char* file_name
  Outputs : no.

  Notes :
      Function call:
        fileio - int split_single_input_file (char* filename).
=============================================================================*/
int check_split_input_file (char *file_name) {
    FILE *fp;
    char line[512];
    int  result = 0;

    /* open file, read first line, check if it is marker for combined file */
    fp = fopen (file_name, "r");
    if (fp == NULL) return -4024;
    if (!feof(fp)) {
        if (fgets (line, 511, fp) != NULL) {  /* Read line. */
            if (strncmp (line, "#<<<BEGIN CONFIGFILE", 20) == 0) {  /* its a combined file! */
#ifdef MPI_PRALLEL
                /* MPI=============================================== */
                if (my_node==ROOT && print_level>=1) printf ("Combined input file detected.\n");
                /* MPI=============================================== */
#else
                if (print_level >= 1) printf ("Combined input file detected.\n");
#endif
                result = 1;
            }
        }
    }
    fclose (fp);

    if (result == 1) split_single_input_file (file_name);  /* its combined, split it! */

    return result;
}

/*=============================================================================
  Function Name : combine_files
  Description   : Combine a number of files into one.
                  The number of strings must be provided!
                  The first string must be the input filename;
                  Then we have alternating:
                  - strings, which are printed into the file;
                  - filenames, which are appended to the file,

  Inputs  : int count, ...
  Outputs : no.

  Notes : no.
=============================================================================*/
int combine_files (int count, ...) {
    va_list ap;

    FILE *fp_out;
    FILE *fp_in;
    int  i;
    char *file_name;  /* to store filename stuff */
    char *line;

    file_name = (char*) malloc (MAX_FILENAME_LENGTH * sizeof (char));
    line = (char*) malloc (2048 * sizeof (char));
    if((file_name == NULL) || (line == NULL)) {
        printf ("Error: insufficient memory!\n");
        return -4025;
    }

    va_start (ap, count);
    if (count < 1) {  /* no output file! exit */
        printf ("Error: too few arguments for file combination!\n");
        va_end (ap);
        return -4026;
    }

    strcpy (file_name, va_arg (ap, char*));  /* obtain filename from next argument */
    fp_out = fopen (file_name, "w");  /* open output file */
    if (fp_out == NULL) {
        printf ("Error! Cannot open combined output file %s.\n", file_name);
        return -4027;
    }

    /* create combined output file */
    for (i=1; i<count; i++){
        if ((i % 2) == 1){  /* it's a string to include */
            strcpy (line, va_arg(ap, char*));  /* obtain string from next argument */
            fprintf (fp_out, "%s\n", line);
        }
        else {  /* it's a filename to include */
            strcpy (file_name, va_arg (ap, char*));  /* obtain filename from next argument */
            /* open file, read lines, copy them into output file */
            fp_in = fopen (file_name, "r");
            if (fp_in == NULL) {printf ("Error: Cannot open file %s for reading.\n", file_name); return -4028;}
            while (!feof (fp_in)) {  /* until file ends */
                /* now the line has been read; put into output file */
                if (fgets (line, 2048, fp_in) != NULL) fputs (line, fp_out);
            }
            fclose (fp_in);
        }
    }
    fclose (fp_out);
    va_end (ap);

    return 0;
}
