/******************************************************************************
  Module Name : aivxyz.c
  Module Date : 05/12/2015
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Output the space-positions (x, y, z) of ions (a), interstitials (i)
                and vacancies (v) in the form of .cfg.

  Others :
      Error numbers in this module 13000-13999.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "aivxyz.h"

/*=============================================================================
  Function Name : init_aivxyz
  Description   : Init aiv_pos.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
int init_aivxyz (void) {

    aiv_i = 0;
    aiv_pos = (struct AIV *) calloc (aiv_num, sizeof (struct AIV));
    if (aiv_pos == NULL) return -13000;  /* cannot allocate memory */

    return 0;
}

/*=============================================================================
  Function Name : aivxyz_output
  Description   : Store A-I-V positions in the form of cfg:
                  x, y, z, type, tab, mater, element.

  Inputs  :
            char *base_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int store_aivxyz_cfg (char *base_name) {
	  FILE *fp;
    int  base_name_length;
    char *str_temp;
    int i;

    base_name_length = strlen (base_name);
    str_temp  = (char*) malloc (sizeof (char) * (base_name_length + 100));
    strncpy (str_temp, base_name, base_name_length + 1);

    strcat (str_temp, "aiv.xyz.cfg");
    if (print_level >= 2) printf ("Storing aiv_xyz to:      %s\n", str_temp);

    fp = fopen (str_temp, "wt");
    if (fp == NULL) return -13001;

    /* .cfg format */
    fprintf (fp, "Number of particles = %i\n", aiv_i);
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
    fprintf (fp, "entry_count = %i\n", 7);
    fprintf (fp, "auxiliary[0] = type (A-I-V)\n");
    fprintf (fp, "auxiliary[1] = tab (A-ion)\n");
    fprintf (fp, "auxiliary[0] = mater\n");
    fprintf (fp, "auxiliary[1] = element\n");
    fprintf (fp, "%i\n", 4);
    fprintf (fp, "aiv_pos\n");

    for (i=0; i<aiv_i; i++)
        fprintf (fp, "%f\t%f\t%f\t%i\t%i\t%i\t%i\n",
        	  aiv_pos[i].x / target_size_x, aiv_pos[i].y / target_size_y, aiv_pos[i].z / target_size_z,
        	  aiv_pos[i].type, aiv_pos[i].tab, aiv_pos[i].mater, aiv_pos[i].element);

    free (aiv_pos);
    fclose (fp);

    return 0;
}
