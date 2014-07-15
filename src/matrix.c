/******************************************************************************
  Module Name : matrix.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Scattering matrix: sin^2(theta_CM/2).

  Others :
      Refers to iradina and Corteo.
      Corteo: Adapted from corteomatrix.c.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "matrix.h"

/* Scattering matrix: sin^2(theta_CM/2) */
float matrix[DIME*DIMS];

/* matrix elements calulation parameters */
#define NSUM  1000  /* number of terms in Gauss-Mehler quadrature sum when computing the matrix */
#define NSUM2 100   /* number of terms in Gauss-Mehler quadrature sum when evaluation the screened
                       cross-section */
#define THETAERR -1000.0

/* parameters used to compute matrix, stored in matrix file */
#define HEADERSIZE 7
float headerRef[HEADERSIZE] = {MINE, DIME, MAXE, MINS, DIMS, MAXS, NSUM};

/* pointer to screening function */
double (*PHI) (double x);

/*=============================================================================
  Function Name : phi_none
  Description   : No screening (screening = 1).

  Inputs  : double x
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double phi_none (double x) {
    return 1.0;
}

/*=============================================================================
  Function Name : phi_univ
  Description   : Universal screening function.

  Inputs  : double x
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double phi_univ (double x) {
    return 0.1818 * exp (-3. * x) + 0.5099 * exp (-0.9423 * x) +
           0.2802 * exp (-0.4028 * x) + 0.02817 * exp (-0.2016 * x);
}

/*=============================================================================
  Function Name : phi_LJ
  Description   : Lenz-Jensen screening function.

  Inputs  : double x
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double phi_LJ (double x) {
    double y = 3.108 * sqrt(x);

    return exp (-y) * (1 + y + 0.3344 * y * y + 0.0485 * y * y * y + 2.647e-3 * y * y * y * y);
}

/*=============================================================================
  Function Name : H
  Description   : Following functions used to compute scattering angle following
                  Gauss-Mehler quadrature.

  Inputs  :
            double u
            double x0
            double epsilon
            double s
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double H (double u, double x0, double epsilon, double s) {
    double x0u = x0 / u;
    double eps1 = 1. / epsilon;
    double ss = s * s;

    return sqrt ((1 - u * u) / (1 - PHI (x0u) / x0u * eps1 - ss / x0u / x0u));
}

/*=============================================================================
  Function Name : func_x0
  Description   : The root of this function (vs x0) is the minimal approach
                  distance at energy epsilon and impact parameter s.

  Inputs  :
            double x0
            double epsilon
            double s
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double func_x0 (double x0, double epsilon, double s) {
    double sx0 = s / x0;

    return 1. - PHI (x0) / (x0 * epsilon) - sx0 * sx0;
}

/*=============================================================================
  Function Name : find_x0
  Description   : Use bisection (Newton) method to find the minimal approach
                  distance x0 (return THETAERR on error).

  Inputs  :
            double epsilon
            double s
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double find_x0 (double epsilon, double s) {
    double func_x0_val1, func_x0_val2, func_x0_valm;

    double temp = 1.0 / (2.0 * epsilon);
    /* inital guesses: Mendenhall & Weller NIMB 58(1991)11, eq. 15 */
    double x2 = temp + sqrt (temp + s * s);
    double x1 = x2 / 10.;

    func_x0_val1 = func_x0 (x1, epsilon, s);  /* should be always negative */
    func_x0_val2 = func_x0 (x2, epsilon, s);  /* should be always positive (starting ~1.0) */
    func_x0_valm = func_x0 ((x1 + x2) / 2., epsilon, s);

    if (func_x0_val1 >= 0.) {
         /* initial guess for x1 too optimistic, start with a safe value (much longer) */
         x1 = MINS;
         func_x0_val1 = func_x0 (x1, epsilon, s);  /* should be always negative */
    }

    if (func_x0_val1>=0. || func_x0_val2<=0.)
        return THETAERR;  /* error, values should be on each side of 0 */

    do {
        if (func_x0_valm < 0.) {
            func_x0_val1 = func_x0_valm;
            x1 = (x1 + x2) / 2.;
        }
        else {
            func_x0_val2 = func_x0_valm;
            x2 = (x1 + x2) / 2.;
        }
        func_x0_valm = func_x0 ((x1 + x2) / 2., epsilon, s);
    } while (fabs (func_x0_valm) > 1e-14);
    /* 1.-XXX wont be more precise than 2^(-52)=2e-16
       but using 1e-10 saves time and is enough precise for float conversion later on */

    return (x1 + x2) / 2.;
}

/*=============================================================================
  Function Name : finds
  Description   : Use bisection method to find the reduced impact parameter s
                  given epsilon and thetaCM (return THETAERR on error).

  Inputs  :
            double epsilon
            double thetaCM
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double finds (double epsilon, double thetaCM) {
    double funcX0val1, funcX0val2, funcX0valm;

    /* inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25 */
    double gamma = (PI - thetaCM) / PI;
    double x0 = find_x0 ((1.0 - gamma * gamma) * epsilon, MINS);
    double x1 = 0.7 * gamma * x0;
    double x2 = 1.0 / (2.0 * epsilon * tan (thetaCM / 2.0));

    funcX0val1 = thetaCM - cal_theta (epsilon, x1, NSUM2);  /* should be always negative */
    funcX0val2 = thetaCM - cal_theta (epsilon, x2, NSUM2);  /* should be always positive (starting ~1.0) */
    funcX0valm = thetaCM - cal_theta (epsilon, (x1 + x2) / 2., NSUM2);

    if (funcX0val1>=0. || funcX0val2<=0.)
        return THETAERR;  /* error, values should be on each side of 0 */

    do {
        if (funcX0valm < 0.) {
            funcX0val1 = funcX0valm;
            x1 = (x1 + x2) / 2.;
        }
        else {
            funcX0val2 = funcX0valm;
            x2 = (x1 + x2) / 2.;
        }
        funcX0valm = thetaCM - cal_theta (epsilon, (x1 + x2) / 2., NSUM2);
    } while (fabs (funcX0valm) > 1e-5);

    return (x1 + x2) / 2.;
}

/*=============================================================================
  Function Name : cal_theta
  Description   : Scattering angle calculated using Gauss-Mehler quadrature
                  (return THETAERR on error).

  Inputs  :
            double epsilon
            double s
            unsigned int nsum
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
double cal_theta (double epsilon, double s, unsigned int nsum) {
  double temp, sum = 0.;
  double x0 = find_x0 (epsilon, s);
  unsigned int j;
  double theta;

  if (x0 == THETAERR) return THETAERR;

  for (j=0; j<nsum/2; j++) {
      temp = H (cos ((2. * j - 1.) * PI / (2. * nsum)), x0, epsilon, s);
      sum += temp;
  }
  theta = 2. * PI * s / nsum / x0 * sum;

  return PI - theta;
}

/*=============================================================================
  Function Name : calc_matrix
  Description   : Compute all the elements of matrix and write matrix to file
                  'corteo.*.mat'; user sets show_progress!=0 to display the
                  progress of this (long) calculation to the console.
                  return 1 if successful, 0 if not able to write file.

  Inputs  :
            unsigned int screening_type
            int show_progress
  Outputs : no.

  Notes :
      Corteo: Adapted from corteomatrix.c.
=============================================================================*/
int calc_matrix (unsigned int screening_type, int show_progress, char *file_name) {
    unsigned long i, j, n_theta_err = 0;
    double theta, sin_theta_by_2;
    FILE *ofp;

    switch (screening_type) {
    case 0 :
        PHI = phi_none;  /* unscreened potential (screening = 1) */
        break;
    case 1 :
        PHI = phi_univ;  /* Universal screening function */
        break;
    case 2 :
        PHI = phi_LJ;    /* Lenz-Jensen screening function */
        break;
    default:
        fprintf(stderr, "Error, screening function of type %u unknown. Stopping here.",
                screening_type);
        exit(1);
    }

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && show_progress)    /* x% */
        fprintf(stdout, "Computing scattering matrix to be stored in %s\n", file_name);
    /* MPI=============================================== */
#else
    //if (show_progress) fprintf (stdout, "Computing scattering matrix");    / * ... */
    if (show_progress)    /* x% */
        fprintf(stdout, "Computing scattering matrix to be stored in %s\n", file_name);
#endif

    /* compute matrix for each reduced energy, reduced impact parameter pair */
    for (j=0; j<DIME; j++) {
#ifdef MPI_PRALLEL
      /* MPI=============================================== */
        if (my_node==ROOT && show_progress && j%(DIME/10)==0) {
            //fprintf (stdout, ".");  / * ... */
            fprintf(stdout, "\b\b\b\b\b\b%.1f%%", j*100./DIME);  /* x% */
            fflush (stdout);
        }
      /* MPI=============================================== */
#else
        if (show_progress && j%(DIME/10)==0) {
            //fprintf (stdout, ".");  / * ... */
            fprintf(stdout, "\b\b\b\b\b\b%.1f%%", j*100./DIME);  /* x% */
            fflush (stdout);
        }
#endif

        for (i=0; i<DIMS; i++) {
            /* calculations (summation) made using double to decrease numerical noise */
            theta = cal_theta (E_val(j), S_val(i), NSUM);
            if (theta == THETAERR) n_theta_err++;
            sin_theta_by_2 = sin (theta / 2.);

            /* store in matrix sin^2(theta_CM/2) as float */
            set_matrix (i + j * DIMS, (float) (sin_theta_by_2 * sin_theta_by_2));
        }
    }

    /* open file corteo.mat and write matrix, including a header indicating
       the parameters used to compute matrix */
    ofp = fopen (file_name, "wb");
    if (ofp == NULL) return 0;
    fwrite ((void*) headerRef, HEADERSIZE, sizeof(float), ofp);
    fwrite ((void*) matrix,    DIME*DIMS,  sizeof(float), ofp);
    fclose (ofp);

#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node==ROOT && show_progress) {
        fprintf (stdout, " done.\n");
        fflush  (stdout);
    }
    /* MPI=============================================== */
#else
    if (show_progress) {
        fprintf (stdout, " done.\n");
        fflush  (stdout);
    }
#endif
    if (n_theta_err) fprintf (stderr, "%lu error(s) evaluating theta\n", n_theta_err);

    return 1;
}

/*=============================================================================
  Function Name : load_matrix
  Description   : Read scattering matrix from file corteo.mat;
                  return 1 if file read correctly, 0 otherwise.

  Inputs  : void
  Outputs : no.

  Notes : no.
=============================================================================*/
int load_matrix (char file_name[]) {
    unsigned int i;
    FILE *ifp;

    float header[HEADERSIZE];
    unsigned int tail;

    ifp = fopen (file_name, "rb");
    if (ifp == NULL) {
        fprintf(stderr, "Can't open %s\n!", file_name);
        return 0;
    }

    /* file exists, get header */
    fread ((void*) header, HEADERSIZE, sizeof (float), ifp);

    for (i=0; i<HEADERSIZE; i++)
        if (header[i] != headerRef[i]) {
            /* "corteo.mat" header does not correspond to currently defined parameters */
            fclose (ifp);
            fprintf (stderr, "Header element %i does not correspond to current parameters in corteo.mat.\n", i);
            return 0;
        }

    /* read file into matrix */
    fread ((void*) matrix, DIME * DIMS, sizeof (float), ifp);

    /* read control data */
    fread ((void*) (&tail), 1, sizeof (int), ifp);
    fclose (ifp);

    return 1;
}

/*=============================================================================
  Function Name : matrix_i
  Description   : Return element i of matrix.

  Inputs  : unsigned long i
  Outputs : no.

  Notes : no.
=============================================================================*/
float matrix_i (unsigned long i) {
    return matrix[i];
}

/*=============================================================================
  Function Name : set_matrix
  Description   : Set element i of matrix.

  Inputs  :
            unsigned long i
            float val
  Outputs : no.

  Notes : no.
=============================================================================*/
void set_matrix (unsigned long i, float val) {
    matrix[i] = val;

    return;
}

/*=============================================================================
  Function Name : prepare_scattering_matrix
  Description   : Creates the matrix with scattering results and stores it to
                  the structure pointed to by scat_matrix.

  Inputs  :
            struct scattering_matrix* scat_matrix
            float proj_M
            int proj_Z
            float target_M
            int target_Z
  Outputs : no.

  Notes :
      Function adapted from corteo.c.

      function call:
        fromcorteo - float d2f (double val);
                     float sqrtdf (double val);
                     void fill_cos_sin_table (float * cos_table,
                                              float * sin_table, float mr).
=============================================================================*/
int prepare_scattering_matrix (struct scattering_matrix *scat_matrix, float proj_M,
                               int proj_Z, float target_M, int target_Z) {

    scat_matrix->screening_length = d2f (0.1f * SCREENCONST / (pow (proj_Z, 0.23)
                                  + pow (target_Z, 0.23)));  /* in nm! */
    scat_matrix->inv_screening_length = d2f (1.0f / scat_matrix->screening_length);
    scat_matrix->mass_ratio = proj_M / target_M;
    scat_matrix->sqrt_mass_ratio = sqrtdf (proj_M / target_M);
    scat_matrix->kfactor_m = 4.0f * scat_matrix->mass_ratio / ((1.0f + scat_matrix->mass_ratio)
                           * (1.0f + scat_matrix->mass_ratio));  /* k factor without the angle part */
    scat_matrix->red_E_conv = 10.0f * scat_matrix->screening_length / ((1.0f + scat_matrix->mass_ratio)
         * proj_Z * target_Z * E2);  /* we need screening length in Angstroms since E2 is given in eVA */

    /* The maximum reduced impact parameter depends on density and flight length. This
       means that we cannot calculate it as a function of projectile and target element only.
       Thus we need to calculate it somewhere else (as a material property of target). This
       is different from corteo, where the scattering matrices are stored for each target
       layer with known density etc. */

    scat_matrix->cos_scat = (float*) malloc (DIME * DIMS * sizeof (float));
    if ((scat_matrix->cos_scat) == NULL) return -1;  /* cannot allocate memory */
    scat_matrix->sin_scat = (float*) malloc (DIME * DIMS * sizeof (float));
    if ((scat_matrix->sin_scat) == NULL) return -1;  /* cannot allocate memory */

    fill_cos_sin_table (scat_matrix->cos_scat, scat_matrix->sin_scat, scat_matrix->mass_ratio);

    return 0;
}

/*=============================================================================
  Function Name : fill_cos_sin_table
  Description   : Fill the table of cos and sin of the scattering angle in the
                  lab frame given a mass ratio.

  Inputs  :
            float * cos_table
            float * sin_table
            float mr
  Outputs :
            float * cos_table
            float * sin_table

  Notes :
      Corteo: Function adapted from corteo.c.

      function call:
        fromcorteo - float matrix_i (unsigned long i);
                     float d2f (double val).
=============================================================================*/
void fill_cos_sin_table (float *cos_table, float *sin_table, float mr) {
    unsigned long i;
    double sin2_thetaby2, cos_theta, cos_theta_lab, sin_theta_lab;

    /* compute scattering angle components */
    for (i=0; i<DIME*DIMS; i++) {
        sin2_thetaby2 = matrix_i (i);
        cos_theta = 1. - 2. * sin2_thetaby2;

        if (cos_theta==-1.0 && mr==1.0) {
            cos_theta_lab = 0.0;  /* peculiar case of head-on collision of identical masses */
            sin_theta_lab = 1.0;
        }
        else {
            cos_theta_lab = (cos_theta + mr) / sqrt (1. + 2. * mr * cos_theta + mr * mr);
            sin_theta_lab = sqrt (1. - cos_theta_lab * cos_theta_lab);
        }
        cos_table[i] = d2f (cos_theta_lab);
        sin_table[i] = d2f (sin_theta_lab);

        /* Modefication from Corteo for iran3d, C. Borschel 2011: */
        /* In some rare cases, when cos=1, then sin becomes "Not a Number".
           To prevent this, I will set the sine to 0 in those cases. */
        if (!((sin_table[i]>=-1) && (sin_table[i]<=1))) sin_table[i] = 0;
    }

    return;
}
