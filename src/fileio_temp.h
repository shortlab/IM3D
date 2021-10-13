/******************************************************************************
  Module Name : fileio.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Contains the some helping function for file access, etc.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
      04/01/14            add subroutine store_results_cfg;
******************************************************************************/
#ifndef FILEIO_H
#define FILEIO_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "const.h"
#include "target.h"
#include "im3d.h"
#include "transport.h"
#include "random.h"
#include "utils.h"
#include "material.h"

/* or ------------------------------------*/
//#include <stdarg.h>

//#include "im3d.h"
//#include "transport.h"

/*--------------------------Global variables-------------------------*/
int single_input_file;  /* if 0, then multiple files (normal).
                           if 1, then a single input file is used */

/*-----------------------------Functions-----------------------------*/
/* A general reader for the various input config files. It parses an ini-like file.
   Whenever it finds a new DataBlock it calls the function pointed to by
   *DataBlockReader with the parameter BlockName. When it finds data it calls the
   function pointed to by DataReader with the parameter-name and its value. */
int read_init_file (int (*read_data_block) (char* block_name),
                    int (*read_data) (char* par_name, char* par_value), char* file_name);

/* writes designated array into file */
int write_int_array_to_file (char *file_name, int *source_array, int count, int file_type);
/* writes designated array into file */
int write_double_array_to_file (char *file_name, double *source_array, int count, int file_type);

/* reads a data block from the configuration file */
int read_config_file_data_block (char *block_name);
/* reads a data block from the configuration file */
int read_config_file_data (char *par_name, char *par_value);

/* store the results of the simulation (arrays with distribution of
   implanted ions, defects etc.) */
int store_results_iradina (char *base_name);

/* Interstitials and so on are stored for each element from each material
   separately, but may also be interesting in sum. So this function does
   all the summing up. */
void sum_up_material_arrays (void);

/* gives the depth distributions */
void do_depth_dist_statistics (void);

/* gives the radial distributions */
void do_radial_dist_statistics (void);

/* store array of depth distribution functions */
int store_depth_dist_array (char *file_name);

/* store array of radial distribution functions */
int store_radial_dist_array (char *file_name);

/* store array of transmitted ions */
int store_transmission_array (char *file_name, struct transmitted_ion *trans_array,
                              int tr_pointer);

/* opens file basename+extension, keeps it open */
FILE *open_file_continuous (char *base_name, char *extension);

/* reads a block of Count float values from file FileName starting at Offset
   and puts these in Array */
int read_float_block (char *file_name, int offset, int count, float *array);

/* does what you think it does */
int write_string_to_file (char *file_name, char *str);

/* prints the contents of a text file to the std out */
int display_a_file (char *file_name);

/* if a single input file is provided, split it */
int split_single_input_file (char *file_name);

/* checks if the given config file is a single combined input file (incl. structure, mat, etc...).
   If so, the file is split into four separate temp files for further "conventional" processing */
int check_split_input_file (char *file_name);

/* combine a number of files into one. */
int combine_files (int count, ...);

#endif
