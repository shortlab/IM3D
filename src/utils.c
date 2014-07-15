/******************************************************************************
  Module Name : utils.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Contains some utility functions.

  Others :
      Error numbers in this file: 4000-4999.
      Refers to iradina and Corteo.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "utils.h"

/*=============================================================================
  Function Name : handle_cmd_line_options
  Description   : Handles all command line arguments. Returns a value > 0 in case
                  the main program should not perform a simultion but something else:
                  1, 2 : nothing;
                     3 : conversion of material file to element files.

  Inputs  :
            int argc - number of commonds;
            char* argv[] - commonds.
  Outputs : no.

  Notes :
      Function call:
        utils  - int print_help_text();
        fileio - int display_a_file (char* Filename).
=============================================================================*/
int handle_cmd_line_options (int argc, char *argv[]) {
    int i = 0;
    int result = 0;

    while ((++i) < argc) {  /* go through command line arguments */
        if (strcmp (argv[i], "-h") == 0) {
#ifdef MPI_PRALLEL
            /* MPI=============================================== */
            if (my_node == ROOT) print_help_text ();
            /* MPI=============================================== */
#else
            print_help_text ();
#endif
            return 1;  /* so the main function knows not to proceed */
        }

        if (strcmp (argv[i], "-l") == 0) {
            if (display_a_file ("license.txt") != 0)
	            printf ("Error: cannot open license.txt\nSee: http://www.gnu.org/licenses/gpl.html\n.");

            return 2;  /* so the main function knows not to proceed */
        }

        /* option for alternative config file name */
        if (strcmp (argv[i], "-c") == 0) {
            i++;
            if ((i<argc) && (strlen (argv[i])<=1023)) {
                strcpy (ConfigFileName, argv[i]);
            }
            else {
	            printf ("Error: the -c option requires a config file name (of less than 1024 characters)\n");
	            return -4001;
            }
        }

        /* option for setting print level */
        if (strcmp (argv[i],"-p") == 0) {
            i++;
            if (i < argc) {
                if (sscanf (argv[i], "%i", &print_level) != 1) {
	                printf ("Error: cannot read level of verbosity from command line option -p\n");
	                return -4002;
	            }
            }
            else {
	            printf ("Error: the -p option requires a following integer number!\n");
	            return -4003;
            }
        }

        /* option for overriding maximum ion number */
        if (strcmp (argv[i], "-n") == 0) {
            i++;
            if (i < argc) {
                if (sscanf (argv[i], "%i", &override_max_ions) !=1 ) {
                    printf ("Error: cannot read level maximum ion number from command line option -n\n");
                    return -4004;
                }
            }
            else {
                printf ("Error: the -n option requires a following integer number!\n");
                return -4005;
            }
        }

        /* option to start iran3d from another prog */
        if (strcmp (argv[i], "-g") == 0) {
            i++;
            if (i < argc) {
                create_status_file = 1;
                start_id_string = (char*) malloc (129 * sizeof (char));
                strncpy (start_id_string, argv[i], 129);
                start_id_string[129] = '\0';
            }
            else {
                printf ("Error: the -g option requires a following string!\n");
                return -4006;
            }
        }

        /* option to override ion energy */
        if (strcmp (argv[i], "-E") == 0) {
            i++;
            if (i < argc) {
                if (sscanf (argv[i], "%f", &override_energy) != 1) {
                    printf ("Error: cannot read override energy from option -E\n");
                    return -4007;
                }
            }
            else {
                printf ("Error: the -E option requires a following float number!\n");
                return -4008;
            }
        }

        /* do not end program before pressing enter.
           useful to check the memory usage of the program */
        if (strcmp (argv[i], "-w") == 0) {
            wait_before_end = 1;
        }

        /* do not simulate, just estimate memory usage of the program */
        if (strcmp (argv[i],"-m") == 0) {
            mem_usage_only = 1;
        }

        /* print memory usage details */
        if (strcmp (argv[i],"-d") == 0) {
            mem_usage_details = 1;
        }

        /* option for converting a material based definition to element based */
        if (strcmp (argv[i], "-conv") == 0) {
            i++;
            if ((i<argc) && (strlen(argv[i])<=1023)) {
                strcpy (ConversionFileName, argv[i]);
                result = 3;
            }
            else {
                printf ("Error: the -conv option requires an output file name (of less than 1024 characters),\n");
                printf ("where the element definition will be stored.\n");
                return -4009;
            }
        }

        /* creates separate elements for each material */
        if (strcmp (argv[i], "-convsep") == 0) {
            conv_create_separate_elements = 1;
        }
    }

    return result;
}

/*=============================================================================
  Function Name : print_hlep_text
  Description   : Prints the help.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
int print_help_text (void) {
    printf ("Usage: iran3d [OPTIONS]\n");
    printf (" Available options: \n");
    printf (" -h            print this help\n");
    printf (" -l            display license\n");
    printf (" -c FILENAME   specify name of config file. Default: Config.in\n");
    printf (" -p NUMBER     specify how much info to print to console. > 0 means much,\n");
    printf ("               < 0 means little\n");
    printf (" -n NUMBER     sets the maximum number of ions to be simulated to NUMBER.\n");
    printf ("               This option overrides the setting from the config file. \n");
    printf (" -E NUMBER     sets the energy of the ions.\n");
    printf ("               This option overrides the setting from the config file. \n");
    printf (" -w            wait for return key before exiting \n");
    printf (" -m            do not simulate, only estimate memory usage (roughly) \n");
    printf (" -d            print details for memory usage (only useful with -m option) \n");
    printf (" -g ID         generate status file while running \n");
    /*printf(" -conv FILE    Converts material based input files to element based input\n");
      printf("               file. If the input file is combined, then the output file\n");
      printf("               will also be combined and written to FILE. Otherwise, the\n");
      printf("               elements are stored in FILE and '.e' is appended to the\n");
      printf("               new compositon file name.\n");
      printf(" -convsep      Should only be used with conv. Create separate elements for\n");
      printf("               each material.\n"); */

    return 0;
}

/*=============================================================================
  Function Name : make_int_array
  Description   : Read comma-seprated values from string and put them into the
                  int array, which has #count entries.
                  The int array must exist already.

  Inputs  :
            char* values
            int count
  Outputs : int* i_array

  Notes : no.
=============================================================================*/
int make_int_array (char *values, int count, int *i_array) {
    int  i = 1;
    char *temp;

    temp = (char*) malloc (sizeof (char)*32);
    temp = strtok (values, ",");
    if (strlen(values) < 1) return -4010;  /* for safety reasons */
    sscanf (temp, "%i", &(i_array[0]));
    while (((temp=strtok(NULL,","))!=NULL) && (i<count)){
        sscanf (temp, "%i", &(i_array[i]));
        i++;
    }

    return 0;
}

/*=============================================================================
  Function Name : make_float_array
  Description   : Read comma-seprated values from string and put them into the
                  float array, which has #count entries.
                  The float array must exist already.

  Inputs  :
            char* values
            int count
  Outputs : int* f_array

  Notes : no.
=============================================================================*/
int make_float_array (char *values, int count, float *f_array) {
    int  i = 1;
    char *temp;

    temp = (char*) malloc (sizeof (char)*32);
    temp = strtok (values, ",");
    sscanf (temp, "%g", &(f_array[0]));
    while (((temp=strtok(NULL,","))!=NULL) && (i<count)) {
        sscanf (temp, "%g", &(f_array[i]));
        i++;
    }

    return 0;
}

/*=============================================================================
  Function Name : fill_int_zero
  Description   : Fill an int arrays with zeros.

  Inputs  :
            int* array
            int count
  Outputs : int* array

  Notes :
      TODO: can be delete ?
=============================================================================*/
void fill_int_zero (int *array, int count) {
    int i;

    for (i=0; i<count; i++) array[i] = 0;

    return;
}

/*=============================================================================
  Function Name : fill_double_zero
  Description   : Fill a float arrays with zeros.

  Inputs  :
            double *array
            int count
  Outputs : int* array

  Notes :
            TODO: can be delete ?
=============================================================================*/
void fill_double_zero (double *array, int count) {
    int i;

    for (i=0; i<count; i++) array[i] = 0;

    return;
}

/*=============================================================================
  Function Name : add_int_array
  Description   : Adds array source to array dest, both arrays must have
                  count entries.

  Inputs  :
            int* dest
            int* source
            int count
  Outputs : int* dest

  Notes : no.
=============================================================================*/
void add_int_array (int *dest, int *source, int count) {
    int i;

    for (i=0; i<count; i++) dest[i] += source[i];

    return;
}

/*=============================================================================
  Function Name : calculate_normalization_factor
  Description   : For converting units to 1/cm^3 per 1/cm^2.

  Inputs  : int num_of_ions
  Outputs : no.

  Notes : no.
=============================================================================*/
void calculate_normalization_factor (int num_of_ions) {

    if (normalize_output == 1) {  /* x -> z */
        unit_conversion_factor = ((double) (1.0 / (cell_size_x * cell_size_y * cell_size_z * 1e-21)))
                    / (((double) (num_of_ions)) / (ion_vx * target_size_y * target_size_z * 1e-14));
    }
    else {
        unit_conversion_factor = 1.0;
    }

    return;
}

/*=============================================================================
  Function Name : get_leaving_direction
  Description   : Get ions/target atoms leaving direction.

  Inputs  : float vx, flaot vy, float vz
  Outputs : leaving_direction

  Notes : no.
=============================================================================*/
int get_leaving_direction (double vx, double vy, double vz) {
    /* This is necessary. In some extremely rare cases, leaving direction is
       undefined otherwise --> crash! */
    int leaving_direction = 0;

    if (vx >= 0 && vy >= 0 && vz >= 0) leaving_direction = 0;  /* thin film */
    if (vx >= 0 && vy <  0 && vz >= 0) leaving_direction = 1;
    if (vx >= 0 && vy <  0 && vz <  0) leaving_direction = 2;
    if (vx >= 0 && vy >= 0 && vz <  0) leaving_direction = 3;
    if (vx <  0 && vy >= 0 && vz >= 0) leaving_direction = 4;
    if (vx <  0 && vy <  0 && vz >= 0) leaving_direction = 5;
    if (vx <  0 && vy <  0 && vz <  0) leaving_direction = 6;
    if (vx <  0 && vy >= 0 && vz <  0) leaving_direction = 7;

    return leaving_direction;
}

/*=============================================================================
  Function Name : write_status_file
  Description   : Create a file that hols status information on iran3D.
                  Can be used to monitor iran3d's status from another program.

  Inputs  :
            char* status_text
            int ion_number
  Outputs : no.

  Notes : no.
=============================================================================*/
int write_status_file (char *status_text, int ion_number) {
    FILE *f_pointer;

    f_pointer = fopen ("ir_state.dat", "w");
    if (f_pointer == NULL) {
        return -4011;
    }
    else {
        fprintf (f_pointer, "iran3d %i.%i.%i\n", VERSION, SUBVERSION, SUBSUBVERSION);
        fprintf (f_pointer, "%s\n", start_id_string);
        fprintf (f_pointer, "%s\n", status_text);
        fprintf (f_pointer, "%i\n", ion_number);
        fclose  (f_pointer);
        return 0;
    }
}

/*=============================================================================
  Function Name : ZBL_and_deri
  Description   : Returns the largest float that is smaller than the fltInput
                  the following code presumes a IEEE754-conform bit-representation
                  of floats it further presumes that int and float have exactly
                  the same bit-length!

  Inputs  : loat* flt_input
  Outputs : float* flt_output

  Notes : no.
=============================================================================*/
void get_float_one_bit_smaller (float *flt_input, float *flt_output) {
    int temp;

    temp = (*((int*) (flt_input))) - 1;
    (*flt_output) = *((float*) (&(temp)));
}

/*=============================================================================
  Function Name : ignoreline
  Description   : Skip to end of line, including end of line character.

  Inputs  : no.
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoutil.c.
=============================================================================*/
void ignore_line (FILE *ifp) {
    fscanf (ifp, "%*[^\n]");
    fscanf (ifp, "%*1[\n]");
}

/*=============================================================================
  Function Name : d2f
  Description   : Convert double to float while preventing underflow.

  Inputs  : double val
  Outputs : double val

  Notes :
      Corteo: Adapted from corteoutil.c.
=============================================================================*/
float d2f (double val) {
    if (val < 0.0) return -val<FLT_MIN?0.0f:(float)val;

    return val<FLT_MIN?0.0f:(float)val;
}

/*=============================================================================
  Function Name : sqrtdf
  Description   : Return sqrt in float format while preventing underflow.

  Inputs  : no.
  Outputs : double val

  Notes :
      Corteo: Adapted from corteoutil.c.
      Function call:
        fromcorteo - float d2f (double val).
=============================================================================*/
float sqrtdf (double val) {
    return d2f (sqrt (val));
}

/*=============================================================================
  Function Name : a2f
  Description   : Convert a string to float while preventing underflow.

  Inputs  : char * s
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoutil.c.

      Function call:
        fromcorteo - float d2f (double val).
=============================================================================*/
float a2f (char *s) {
    return d2f (atof (s));
}

/* no use now, maybe used in the future */
/*=============================================================================
  Function Name : fast_sqrt
  Description   : Sqrt of val is the product of the sqrt of its mantissa and
                  sqrt of its exponent.

  Inputs  : float val
  Outputs : no.

  Notes :
      WARNING: approximate solution, precise to 0.0015%, use when precision
               is not critical.
=============================================================================*/
float fast_sqrt (float val) {
    if (*(unsigned long *)&val==0) return 0.0f;
    return sqrt_table[((*(unsigned long *)&val)>>7)&0xFFFF]
         * sqrt_table_exp[((*(unsigned long *)&val)>>23)&0xFF];
}

/*=============================================================================
  Function Name : inv_sqrt
  Description   : 1/sqrt of val is the product of the 1/sqrt of its mantissa
                  and 1/sqrt of its exponent.
                  WARNING: approximate solution, precise to 0.0015%, use when
                  precision is not critical.

  Inputs  :
            float val
  Outputs : no.

  Notes : no.
=============================================================================*/
float inv_sqrt (float val) {

    if ((*(unsigned long *)&val) == 0) return 1.0f / val;  // prevent division by 0
    return inv_sqrt_table[((*(unsigned long *)&val)>>7)&0xFFFF]
         * inv_sqrt_table_exp[((*(unsigned long *)&val)>>23)&0xFF];
}

/*=============================================================================
  Function Name : mat_mul
  Description   : Multiply matrix: a[3*3] * b[3*3] = c[3*3].

  Inputs  : double a[][3], double b[][3]
  Outputs : double c[][3]

  Notes : no.
=============================================================================*/
void mat_mul (double a[][3], double b[][3], double c[][3]) {
    int i, j, k;

    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            c[i][j] = 0;
            for (k=0; k<3; k++) c[i][j] += a[i][k] * b[k][j];
        }
    }

    return;
}

/*=============================================================================
  Function Name : mat_mul
  Description   : Multiply matrix: a[3*3] * b[3] = c[3].

  Inputs  : double a[][3], double b[][3]
  Outputs : double c[][3]

  Notes : no.
=============================================================================*/
void mat_mul2 (double a[][3], double b[], double c[]) {
    int i, j;

    for (i=0; i<3; i++) {
        c[i] = 0;
        for (j=0; j<3; j++) {
            c[i] += a[i][j] * b[i];
        }
    }

    return;
}

/*=============================================================================
  Function Name : dot_product
  Description   :

  Inputs  :
            double a[]
            double b[]
  Outputs : double c

  Notes : no.
=============================================================================*/
double dot_product (double a[], double b[]) {

    double c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    return (c);
}

/*=============================================================================
  Function Name : cross_product
  Description   :

  Inputs  :
            double a[]
            double b[]
  Outputs : double c[]

  Notes : no.
=============================================================================*/
void cross_product (double a[], double b[], double c[]) {

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];

    return;
}

/*=============================================================================
  Function Name : copy_int_array
  Description   :

  Inputs  :
            int a[], b[]
            int num_elements
  Outputs : double sum

  Notes : no.
=============================================================================*/
 int copy_int_array (int a[], int b[], int num_elements) {
    int i;

    for (i=0; i<num_elements; i++) b[i] = a[i];

    return 0;
}

/*=============================================================================
  Function Name : copy_double_array
  Description   :

  Inputs  :
            double a[], b[]
            int num_elements
  Outputs : double sum

  Notes : no.
=============================================================================*/
 double copy_double_array (double a[], double b[], int num_elements) {
    int i;

    for (i=0; i<num_elements; i++) b[i] = a[i];

    return 0;
}

/*=============================================================================
  Function Name : sum_array
  Description   :

  Inputs  :
            double a[]
            int num_elements
  Outputs : double sum

  Notes : no.
=============================================================================*/
 double sum_array (double a[], int num_elements) {
    int i;
    double sum = 0;

    for (i=0; i<num_elements; i++) sum = sum + a[i];

    return (sum);
}

/*=============================================================================
  Function Name : max_loc_abs
  Description   :

  Inputs  :
            double a[]
            int n
  Outputs : int ML

  Notes : no.
=============================================================================*/
int max_loc_abs (double a[], int n) {
    int i, ML;
    double temp, max;

    ML = 0;
    max = fabs (a[0]);
    for (i=0; i<n; i++) {
        temp = fabs (a[i]);
        if (max < temp) {
            max = temp;
            ML = i;
        }
    }

    return (ML);
}
