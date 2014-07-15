/***********************************************************************
  Module Name : utils.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Contains some auxiliary functions.

  Others :
      Refers to iradina and Corteo.

      Revision History:
      Date    Rel Ver.    Notes
***********************************************************************/
#ifndef UTILS_H
#define UTILS_H

/*------------------------------Includes-----------------------------*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "fileio.h"
#include "iran3d.h"
#include "transport.h"
#include "target.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include <float.h>

//#include "fileio.h"

/*------------------------------Defines------------------------------*/
#define MIN(x,y) (x<y)?x:y
#define MAX(x,y) (x>y)?x:y

/*--------------------------Global variables-------------------------*/
/* if 1 then, separate elements are created for each material,
   so some elements may appear more than once */
int conv_create_separate_elements;

/*-----------------------------Functions-----------------------------*/
/* ... */
int handle_cmd_line_options (int argc, char *argv[]);

int print_help_text (void);

/* read comma-seprated values from string and put them into
   the int array, which has #count entries */
int make_int_array (char *values, int count, int *i_array);

/* read comma-seprated values from string and put them into
   the float array, which has #count entries. The float array
   must exist already */
int make_float_array (char *values, int count, float *f_array);

/* fills an int array with zeros */
void fill_int_zero (int *array, int count);

/* fills a double array with zeros */
void fill_double_zero (double *array, int count);

/* adds array source to array dest */
void add_int_array (int *dest, int *source, int count);

/* for converting units to 1/cm^2 per 1/cm^2 */
void calculate_normalization_factor (int num_of_ions);

/* get ions/target atoms leaving direction */
int get_leaving_direction (double vx, double vy, double vz);

/* creates a file, that describes iran3d's running state */
int write_status_file (char *status_text, int ion_number);

/* returns the largest float that is smaller than the fltInput */
void get_float_one_bit_smaller (float *flt_input, float *flt_output);

void ignore_line (FILE *ifp);
float d2f (double val);
float sqrtdf (double val);
float a2f (char *s);

float inv_sqrt (float val);
float fast_sqrt (float val);

/* multiply matrix: a[m*n] * b[n*k] = c[m*k]. */
//void mat_mul (double a[], double b[], int m, int n, int k, double c[]);
void mat_mul (double a[][3], double b[][3], double c[][3]);
void mat_mul2 (double a[][3], double b[], double c[]);

double dot_product (double a[], double b[]);

void cross_product (double a[], double b [], double c[]);

int copy_int_array (int a[], int b[], int num_elements);

double copy_double_array (double a[], double b[], int num_elements);

double sum_array (double a[], int num_elements);

int max_loc_abs (double a[], int n);

#endif
