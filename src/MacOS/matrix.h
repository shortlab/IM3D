/******************************************************************************
  Module Name : matrix.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Scattering matrix: sin^2(theta_CM/2).

  Others :
      Refers to iradina and Corteo
      Corteo: Adapted from corteomatrix.c.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

/*------------------------------Includes-----------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "index.h"
#include "const.h"
#include "utils.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include <stdlib.h>

//#include "const.h"
//#include "index.h"
//#include "material.h"

/*------------------------------Defines------------------------------*/
#define MAX_ELEMENT_NO 93  /* maximum number of elements in simulation +1 */

/*--------------------------Global variables-------------------------*/
struct scattering_matrix  /* This structure holds information for a scattering events of two
                             particles, for each possible energy and impact parameter. */
{
    float screening_length;      /* screening length from ZBL85, p45, eq.2-60, but in [nm] */
    float inv_screening_length;  /* 1/screening length in 1/nm */
    float mass_ratio;            /* ... */
    float sqrt_mass_ratio;       /* we will need this occasionally */
    float kfactor_m;             /* mass part of the kinematic factor.
				                            This is called EC in the TRIM code. */
    float red_E_conv;            /* reduced energy conversion factor. This is called FI in TRIM */
    float *cos_scat;             /* cosines of scattering angles, for each energy and p */
    float *sin_scat;             /* sines of scattering angles, for each energy and p */
};

/* The following structures hold scattering matrices for each possible projectile and
   target combination.
   The ion gets its own matrix, because the ion might be a specific isotope which is not
   the naturally most abundant.
   For the other elements, we will assume most abundant isotopes only, because memory is
   not unlimited.
   In case we need to care about specific target isotopes, the MAGIC algorithm should be
   used instead of the ScatBase approach. */
struct scattering_matrix ion_scattering_matrix[MAX_ELEMENT_NO];
/* First index: projectile, second: target */
struct scattering_matrix scattering_matrices[MAX_ELEMENT_NO][MAX_ELEMENT_NO];

/*-----------------------------Functions-----------------------------*/
int calc_matrix (unsigned int screening_type, int show_progress, char *file_name);
int load_matrix (char file_name[]);
float matrix_i (unsigned long i);
void set_matrix (unsigned long i, float val);
double cal_theta (double epsilon, double s, unsigned int nsum);

/*Creates the matrix with scattering results and stores it to the structure pointed to by ScatMatrix*/
/* This function needs to be here, to exclude some warning */
int prepare_scattering_matrix (struct scattering_matrix *scat_matrix, float proj_M, int proj_Z,
                               float target_M, int target_Z);

/* corteo function to calc cos and sin of scattering angles, adapted from corteo.h */
void fill_cos_sin_table (float *cos_table, float *sin_table, float mr);

#endif
